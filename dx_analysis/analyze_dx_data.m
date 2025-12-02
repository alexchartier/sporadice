%% analyze_dx_data.m
% Analysis of DX maps data provided by Chris Deacon

%% Set inputs
in_fn_1 = '~/data/sporadic_e/NA Es spots 28 MHz 11 May 2024 0000 - 0400.xlsx';
in_fn_2 = '~/data/sporadic_e/NA Es spots 50 MHz 11 May 2024 0000 - 0400.xlsx';
plt_fn_fmt = '~/superdarn/plots/20240511/esmaps/esmaps_{yyyymmdd_HHMM}.png';

times = datenum(2024, 5, 11):5/60/24:datenum(2024, 5, 11, 4, 0, 0);
lats = 0:2:60;
lons = -140:2:-60;
cornell_box = [44, 47, -81, -76];

Re = 6371E3;
hmE = 100E3;

%% Load
data = [readtable(in_fn_1); readtable(in_fn_2)];

%% Calculate coordinates
for ind = 1:height(data)
    txloc = regexprep(data.DX{ind},{'.*\(','\).*'},'');
    rxloc = regexprep(data.Spotter{ind},{'.*\(','\).*'},'');
    [data.txlat(ind), data.txlon(ind)] = maidenhead_to_ll(txloc);
    [data.rxlat(ind), data.rxlon(ind)] = maidenhead_to_ll(rxloc);
end    

%% Calculate sporadic E min density on each ray
[data.lat_mid, data.lon_mid] = ...
    midpointLatLon(data.txlat, data.txlon, data.rxlat, data.rxlon);
data.gc_dist = greatcircle(data.txlat, data.txlon, data.rxlat, data.rxlon);

[az, elev, rg] = geodetic2aer( ...
    data.txlat, data.txlon, 0, data.lat_mid, data.lon_mid, hmE, wgs84Ellipsoid, 'degrees');
data.interaction_angle = elev;

data.nmE_min = freq2elec(data.QRG .* cosd(elev + 90));

% fix lons
data.txlon(data.txlon > 180) = data.txlon(data.txlon > 180) - 360; 
data.rxlon(data.rxlon > 180) = data.txlon(data.txlon > 180) - 360; 
data.lon_mid(data.lon_mid > 180) = data.lon_mid(data.lon_mid > 180) - 360; 

%% gridwise sporadic E
nmE_grid = zeros(length(times), length(lats), length(lons));
lat_d = lats(2) - lats(1);
lon_d = lons(2) - lons(1);
time_d = times(2) - times(1);
for ind = 1:height(data)
    [lat, latd, lati] = closest(lats, data.lat_mid(ind));
    [lon, lond, loni] = closest(lons, data.lon_mid(ind));
    [time, td, ti] = closest(times, datenum(data.UTCDate_Time(ind)));
    assert(latd  < lat_d, 'too far')
    assert(lond  < lon_d, 'too far')
    assert(td  < time_d, 'too long')
   
    if nmE_grid(ti, lati, loni) < data.nmE_min(ind)
        nmE_grid(ti, lati, loni) = data.nmE_min(ind);
    end
end

nmE_grid(nmE_grid == 0) = NaN;

%% Plotting
for t = 1:length(times)
    hold on
    % setup
    %m_proj('lambert','long',[-80 -66],'lat',[36 50]);
    m_proj('lambert','long',[-120 -66],'lat',[24 55]);
    m_coast('patch',[1 .85 .7]);
    m_grid('box','fancy','tickdir','in');
    
    % NmEs
    m_pcolor(lons - lon_d/2, lats - lat_d/2, squeeze(nmE_grid(t, :, :)))

    % links
    ti = abs(datenum(data{:, 1} - times(t))) < (times(2) - times(1)) / 2;
    txlat = data{ti, 13};
    txlon = data{ti, 14};
    rxlat = data{ti, 15};
    rxlon = data{ti, 16};
    for i = 1:length(txlat)
        [lat, lon, ~] = greatcircle(txlat(i), txlon(i), rxlat(i), rxlon(i));
        m_plot(lon, lat, '-r', 'LineWidth', 2)
    end

    m_plot(...
        [cornell_box(3:4), cornell_box(4), cornell_box(3), cornell_box(3)], ...
        [cornell_box(1), cornell_box(1:2), cornell_box(2), cornell_box(1)], ...
        '-m', 'LineWidth', 3 ...
        )
    title(datestr(times(t), 'YYYY-mm-dd HH:MM'))
    a = colorbar;
    ylabel(a,'NmE (10^{12} el. m^{3})','FontSize',16,'Rotation',90);

    clim([0, 1.5])
    hold off
    pause(0.1)

    export_fig(filename(plt_fn_fmt, times(t)))
    clf
end



%% timeseries plot

boxi = data{:, 17} > cornell_box(1) & data{:, 17} < cornell_box(2) & ...
    data{:, 18} > cornell_box(3) & data{:, 18} < cornell_box(4);
plot(data{boxi, 1}, data{boxi, 21}, '.', 'MarkerSize', 30);
title('Ham VHF link estimates of NmE between 40-47 N, 76-81 W')
ylabel('NmE (10^{12} el. m^{-3})')
ylim([0, 2.2])
xlabel('Time (UT)')
grid on
grid minor
set(gca, 'FontSize', 20)




































