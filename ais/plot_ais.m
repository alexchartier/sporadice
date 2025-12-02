%% Set inputs
delta_t = 10/60/24;
times = datenum(2024, 5, 9):delta_t:datenum(2024, 5, 12);
in_fn_fmt = '~/data/sporadice/ais_proc_l2/{yyyymmdd_HHMM}.mat';
station_fn = '~/data/sporadice/ais_station_info/station_info.mat';
era5_fn_fmt = '~/data/sporadice/era5/max_duct_strength_{yyyymmdd_HHMM}.nc';
plt_fn_fmt = '~/data/sporadice/plots/ais_era5_maps/{yyyymmdd_HHMM}.png';

%% Plot links
station_info = load_station_info(station_fn);
cols = wide_range_colormap(24);
colormap([1 1 1; cols]);
caxis([-30 95]);
era_prev = [];
era_next = [];
era_prev_time = NaN;
era_next_time = NaN;
for t = 1%:length(times)
    this_time = times(t);
    % load each hour
    clf
    if (this_time * 24) == round(this_time * 24)
        ais = loadstruct(filename(in_fn_fmt, this_time));
        
        idx_es = ais.link_classification == 'other';
        idx_tropo = ais.link_classification == 'tropo';
        ais_es = ais(idx_es, :);
        ais_tropo = ais(idx_tropo, :);
    end

    base_time = floor(this_time * 24) / 24;
    next_time = base_time + 1/24;

    % Interpolate ERA5 between hourly files to match the plot time
    if base_time ~= era_prev_time
        era_prev = load_era5_if_exists(era5_fn_fmt, base_time);
        era_prev_time = base_time;
    end

    if next_time ~= era_next_time
        era_next = load_era5_if_exists(era5_fn_fmt, next_time);
        era_next_time = next_time;
    end

    if ~isempty(era_prev) && ~isempty(era_next)
        w = (this_time - base_time) / (next_time - base_time);
        era_lon = era_prev.lon;
        era_lat = era_prev.lat;
        max_duct_strength = (1 - w) .* era_prev.max_duct_strength + ...
            w .* era_next.max_duct_strength;
    elseif ~isempty(era_prev)
        era_lon = era_prev.lon;
        era_lat = era_prev.lat;
        max_duct_strength = era_prev.max_duct_strength;
    elseif ~isempty(era_next)
        era_lon = era_next.lon;
        era_lat = era_next.lat;
        max_duct_strength = era_next.max_duct_strength;
    else
        warning('Missing ERA5 grids for %s', datestr(this_time));
        continue
    end

    tidx_es = abs(ais_es.time - this_time) < (delta_t / 2);
    tidx_tropo = abs(ais_tropo.time - this_time) < (delta_t / 2);
    C = load('coastlines');
    nans_es = nan(size(ais_es.rxlat(tidx_es)));
    nans_tropo = nan(size(ais_tropo.rxlat(tidx_tropo)));
    latarr_es = [ais_es.txlat(tidx_es), ais_es.rxlat(tidx_es), nans_es]';
    lonarr_es = [ais_es.txlon(tidx_es), ais_es.rxlon(tidx_es), nans_es]';
    latarr_tropo = [ais_tropo.txlat(tidx_tropo), ais_tropo.rxlat(tidx_tropo), nans_tropo]';
    lonarr_tropo = [ais_tropo.txlon(tidx_tropo), ais_tropo.rxlon(tidx_tropo), nans_tropo]';
    hold on
    imagesc(era_lon, era_lat, max_duct_strength)
    set(gca, 'ydir', 'normal', 'Layer', 'top');
    hc = colorbar;
    set(hc,'FontSize', 14, 'ytick', -25:5:95);
    hcl = get(hc, 'ylabel');
    set(hcl, 'string', 'M-Excess ( \Delta M)', 'FontSize', 14)
    clim([-25, 95]);
    h3 = plot(station_info.lon, station_info.lat, '.k', 'MarkerSize', 10, ...
        'DisplayName', 'Rx Station');
    h1 = plot(lonarr_es(:), latarr_es(:), '-r', 'LineWidth', 4, ...
        'DisplayName', 'Potential Es');
    h2 = plot(lonarr_tropo(:), latarr_tropo(:), '-k', 'LineWidth', 1, ...
        'DisplayName', 'Tropo');
    plot(C.coastlon, C.coastlat,'k')
    ylabel(['Latitude (' char(176) 'N)'],'FontSize',14);
    xlabel(['Longitude (' char(176) 'E)'],'FontSize',14);
    legend([h2, h1, h3])
    title(datestr(this_time, 'yyyy/mm/dd HH:MM'))
    hold off
    xlim([-180, -50])
    ylim([10, 80])
    grid on
    grid minor
    
    export_fig(filename(plt_fn_fmt, this_time))
    pause(0.1)

end

function era5 = load_era5_if_exists(era5_fn_fmt, this_time)
era_fn = filename(era5_fn_fmt, this_time);
if isfile(era_fn)
    era5 = load_nc(era_fn);
else
    era5 = [];
end
end
