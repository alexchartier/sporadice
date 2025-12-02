%% ais_l1_to_l2.m
% quality flag AIS output using the link data + ERA5 duct strengths

clear
%% Configuration
times = datenum(2024, 5, 9):1/24:datenum(2024, 5, 12);
% times = datenum(2024, 5, 9):1/24:datenum(2024, 5, 9, 1, 0, 0);

in_fn_fmt = '~/data/sporadice/ais_proc_l1/{yyyymmdd_HHMM}.mat';
out_fn_fmt = '~/data/sporadice/ais_proc_l2/{yyyymmdd_HHMM}.mat';
era5_fn_fmt = '~/data/sporadice/era5/max_duct_strength_{yyyymmdd_HHMM}.nc';
station_dist_thresh_km = 100;
reflect_alt_km = 110;
duct_strength_threshold = -5;
duct_radius_km = 100;
ais_freq_MHz = 162;
station_fn = '~/data/sporadice/ais_station_info/station_info.mat';


%% Loop over times
station_info = load_station_info(station_fn);

% Parallelize the per-hour processing to speed up the L1 -> L2 conversion
if isempty(gcp('nocreate'))
    parpool; 
end
parfor t = 1:length(times)
    process_single_time(...
        times(t), in_fn_fmt, out_fn_fmt, era5_fn_fmt, ...
        reflect_alt_km, ais_freq_MHz, duct_radius_km, ...
        station_info, station_dist_thresh_km);
end


%% Classify the links 
for t = 1:length(times)
    ais = loadstruct(filename(out_fn_fmt, times(t)));

    % Classify based on flags
    ais = classify_links(ais, duct_strength_threshold);
    
    % Save
    out_fn = filename(out_fn_fmt, times(t));
    save(out_fn, 'ais', '-v7');
    fprintf('Saved to %s\n', out_fn)
end


%% Helper function to process one time
function process_single_time(this_time, in_fn_fmt, out_fn_fmt, era5_fn_fmt, ...
    reflect_alt_km, ais_freq_MHz, duct_radius_km, station_info, station_dist_thresh_km)
fprintf('%s\n', datestr(this_time))

% Load AIS table
ais = loadstruct(filename(in_fn_fmt, this_time));

if isempty(ais)
    fprintf('    AIS table is empty, skipping %s\n', datestr(this_time));
    return
end

% Geometry to estimate electron density
ais = get_midpoints(ais);
ais = add_midpoint_elevation(ais, reflect_alt_km);
ais = add_implied_nmax(ais, ais_freq_MHz);

% Calculate ERA5 duct strength within a configurable radius (default 500 km)
ais.max_duct_strength = era5_max_duct_within_radius(...
    filename(era5_fn_fmt, this_time), ...
    ais.txlat, ais.txlon, ais.rxlat, ais.rxlon, duct_radius_km);

% Add quality flags
ais = add_pseudo_station_id(ais, station_info);
ais = add_station_proximity_flags(ais, station_info, station_dist_thresh_km);
ais = add_quality_flags(ais);

% Save
out_fn = filename(out_fn_fmt, this_time);
save(out_fn, 'ais', '-v7');
end

%% Helper to get maximum duct strength within a chosen radius of each link
function duct_vals = era5_max_duct_within_radius(...
    era_fn, txlat, txlon, rxlat, rxlon, radius_km)

if nargin < 6 || isempty(radius_km)
    radius_km = 500;
end
radius_km = double(radius_km);
[lat_vec, lon_vec, duct_grid] = load_era5_grid(era_fn);

txlat = txlat(:);
txlon = txlon(:);
rxlat = rxlat(:);
rxlon = rxlon(:);
txlon(txlon > 180) = txlon(txlon > 180) - 360;
rxlon(rxlon > 180) = rxlon(rxlon > 180) - 360;

n_links = numel(txlat);
duct_vals = nan(n_links, 1);
dist_km = greatcircle(txlat, txlon, rxlat, rxlon);

% Coarser sampling for speed; planar distance approximation
step_km = 100;
buffer_km = radius_km;
buffer_lat_deg = buffer_km / 111; % approximate degrees latitude

for i = 1:n_links
    if isnan(txlat(i)) || isnan(txlon(i)) || isnan(rxlat(i)) || isnan(rxlon(i)) || isnan(dist_km(i))
        continue
    end

    n_samples = max(ceil(dist_km(i) / step_km) + 1, 2);
    [lat_path, lon_path] = greatcircle(txlat(i), txlon(i), rxlat(i), rxlon(i), n_samples);
    lon_path = mod(lon_path + 180, 360) - 180; % wrap to [-180, 180]
    lon_path_unw = rad2deg(unwrap(deg2rad(lon_path)));

    mean_lat = mean([txlat(i), rxlat(i)], 'omitnan');
    lon_scale = max(cosd(mean_lat), 0.2);
    buffer_lon_deg = buffer_km / (111 * lon_scale);

    lat_min = max(min(lat_path) - buffer_lat_deg, min(lat_vec));
    lat_max = min(max(lat_path) + buffer_lat_deg, max(lat_vec));
    lon_min = min(lon_path_unw) - buffer_lon_deg;
    lon_max = max(lon_path_unw) + buffer_lon_deg;

    lon_center = (lon_min + lon_max) / 2;
    lon_vec_unw = lon_vec + 360 * round((lon_center - lon_vec) / 360);

    lat_idx = lat_vec >= lat_min & lat_vec <= lat_max;
    lon_idx = lon_vec_unw >= lon_min & lon_vec_unw <= lon_max;
    if ~any(lat_idx) || ~any(lon_idx)
        continue
    end

    sub_grid = duct_grid(lat_idx, lon_idx);
    lat_sub = lat_vec(lat_idx);
    lon_sub = lon_vec_unw(lon_idx);

    [lon_mesh, lat_mesh] = meshgrid(lon_sub, lat_sub);
    lat_vec_flat = lat_mesh(:);
    lon_vec_flat = lon_mesh(:);

    mean_path_lat = mean(lat_path, 'omitnan');
    lon_scale = max(cosd(mean_path_lat), 0.2);

    dlat_km = (lat_path(:) - lat_vec_flat.').*111;
    dlon_km = (lon_path_unw(:) - lon_vec_flat.').*111*lon_scale;
    dist2 = dlat_km.^2 + dlon_km.^2;

    min_dist2 = min(dist2, [], 1);
    close_idx = min_dist2 <= buffer_km^2;
    if ~any(close_idx)
        continue
    end

    candidates = sub_grid(:);
    duct_vals(i) = max(candidates(close_idx));
end
end

%% Helper to load and standardize an ERA5 duct grid
function [lat_vec, lon_vec, duct_grid] = load_era5_grid(era_fn)
lat_vec = double(ncread(era_fn, 'lat'));
lon_vec = double(ncread(era_fn, 'lon'));
duct_grid = double(ncread(era_fn, 'max_duct_strength'));

duct_grid(duct_grid > 1e19) = -45;
% duct_grid(isnan(duct_grid)) = -45;

% NetCDF is lon x lat; transpose to lat x lon if needed
if size(duct_grid, 1) == numel(lon_vec) && size(duct_grid, 2) == numel(lat_vec)
    duct_grid = duct_grid.';
elseif size(duct_grid, 1) ~= numel(lat_vec) || size(duct_grid, 2) ~= numel(lon_vec)
    error('Unexpected ERA5 dimensions in %s', era_fn);
end

% Standardize coordinates
if any(lon_vec > 180)
    lon_vec(lon_vec > 180) = lon_vec(lon_vec > 180) - 360;
end
[lon_vec, lon_order] = sort(lon_vec);
duct_grid = duct_grid(:, lon_order);

if ~issorted(lat_vec)
    [lat_vec, lat_order] = sort(lat_vec);
    duct_grid = duct_grid(lat_order, :);
end
end

%% Helper to assign a pseudo-station ID based on receiver coordinates
function ais = add_pseudo_station_id(ais, station_info)
assert(all(ismember({'rxlat', 'rxlon'}, ais.Properties.VariableNames)), ...
    'Missing rxlat/rxlon in AIS table.');

rx_lat = ais.rxlat;
rx_lon = ais.rxlon;
rx_lon(rx_lon > 180) = rx_lon(rx_lon > 180) - 360;

% Round to mitigate tiny floating differences
rx_lat = round(rx_lat * 1e3) / 1e3;
rx_lon = round(rx_lon * 1e3) / 1e3;

valid = ~isnan(rx_lat) & ~isnan(rx_lon);
pseudo_id = nan(height(ais), 1);
if ~any(valid)
    ais.rx_station_id = pseudo_id;
    return
end

rx_coords = [rx_lat(valid), rx_lon(valid)];
base_coords = [];

if ~isempty(station_info)
    base_lat = station_info.lat;
    base_lon = station_info.lon;
    base_lon(base_lon > 180) = base_lon(base_lon > 180) - 360;
    base_lat = round(base_lat * 1e3) / 1e3;
    base_lon = round(base_lon * 1e3) / 1e3;
    base_coords = [base_lat(:), base_lon(:)];
end

if isempty(base_coords)
    [~, ~, ids] = unique(rx_coords, 'rows', 'stable');
    pseudo_id(valid) = ids;
else
    [~, idx] = ismember(rx_coords, base_coords, 'rows');
    missing = idx == 0;
    if any(missing)
        new_coords = unique(rx_coords(missing, :), 'rows', 'stable');
        base_coords = [base_coords; new_coords];
        [~, idx] = ismember(rx_coords, base_coords, 'rows');
    end
    pseudo_id(valid) = idx;
end

ais.rx_station_id = pseudo_id;
end

%% Helper to mark stations with neighbours within the distance threshold
function ais = add_station_proximity_flags(ais, station_info, dist_thresh_km)
assert(all(ismember({'rxlat', 'rxlon'}, ais.Properties.VariableNames)), ...
    'Missing rxlat/rxlon in AIS table.');

station_ids = nan(height(ais), 1);
if ismember('rx_station_id', ais.Properties.VariableNames)
    station_ids = ais.rx_station_id;
end

rx_lat = ais.rxlat;
rx_lon = ais.rxlon;
rx_lon(rx_lon > 180) = rx_lon(rx_lon > 180) - 360;

% Round to dampen tiny differences
rx_lat = round(rx_lat * 1e3) / 1e3;
rx_lon = round(rx_lon * 1e3) / 1e3;

nrows = height(ais);
near_in_hour = false(nrows, 1);
near_in_station_info = false(nrows, 1);

valid = ~isnan(rx_lat) & ~isnan(rx_lon) & ~isnan(station_ids);
if ~any(valid)
    ais.rx_station_nearby_in_hour = near_in_hour;
    ais.rx_station_nearby_in_station_info = near_in_station_info;
    return
end

% Representative coordinate for each station in the hour
[unique_ids, ~, ~] = unique(station_ids(valid));
num_stations = numel(unique_ids);
station_coords = nan(num_stations, 2);
for i = 1:num_stations
    idx = valid & station_ids == unique_ids(i);
    station_coords(i, 1) = mean(rx_lat(idx));
    station_coords(i, 2) = mean(rx_lon(idx));
end
station_coords = round(station_coords * 1e3) / 1e3;

% Station info coordinates (if available)
info_lat = [];
info_lon = [];
if ~isempty(station_info)
    info_lat = station_info.lat(:);
    info_lon = station_info.lon(:);
    info_lon(info_lon > 180) = info_lon(info_lon > 180) - 360;
    info_lat = round(info_lat * 1e3) / 1e3;
    info_lon = round(info_lon * 1e3) / 1e3;
    valid_info = ~isnan(info_lat) & ~isnan(info_lon);
    info_lat = info_lat(valid_info);
    info_lon = info_lon(valid_info);
end

for i = 1:num_stations
    this_lat = station_coords(i, 1);
    this_lon = station_coords(i, 2);
    if any(isnan([this_lat, this_lon]))
        continue
    end

    % Check other stations in this hourly file
    dist_vec = greatcircle(this_lat * ones(num_stations, 1), ...
        this_lon * ones(num_stations, 1), ...
        station_coords(:, 1), station_coords(:, 2));
    dist_vec(i) = inf; % ignore self
    has_peer = any(dist_vec <= dist_thresh_km);

    % Check against full station_info list
    has_info_neighbor = false;
    if ~isempty(info_lat)
        dist_info = greatcircle(this_lat * ones(numel(info_lat), 1), ...
            this_lon * ones(numel(info_lon), 1), info_lat, info_lon);
        % Ignore self if this station maps to a station_info row
        self_idx = unique_ids(i);
        if self_idx >= 1 && self_idx <= numel(info_lat)
            dist_info(self_idx) = inf;
        end
        has_info_neighbor = any(dist_info <= dist_thresh_km);
    end

    row_mask = station_ids == unique_ids(i);
    near_in_hour(row_mask) = has_peer;
    near_in_station_info(row_mask) = has_info_neighbor;
end

ais.rx_station_nearby_in_hour = near_in_hour;
ais.rx_station_nearby_in_station_info = near_in_station_info;
ais.bad_receiver = not(ais.rx_station_nearby_in_hour) & ...
    ais.rx_station_nearby_in_station_info;
end

%% Helper to add quality/cleanliness flags
function ais = add_quality_flags(ais)
nrows = height(ais);
var_names = ais.Properties.VariableNames;

ais.bad_position = false(nrows, 1);
ais.aircraft = false(nrows, 1);
ais.likely_corrupt_lon = false(nrows, 1);
ais.known_bad_ship = false(nrows, 1);
ais.bogus_mmsi = false(nrows, 1);

has_mmsi = ismember('mmsi', var_names);
has_tx = all(ismember({'txlat', 'txlon'}, var_names));
has_rx = all(ismember({'rxlat', 'rxlon'}, var_names));

% Flag zeroed positions
if has_tx
    ais.bad_position = (ais.txlat == 0) | (ais.txlon == 0);
end

if has_mmsi
    mmsi = ais.mmsi;
    ais.aircraft = mmsi > 1e8 & mmsi < 2e8;
    ais.known_bad_ship = mmsi == 311001094;

    bogus_list = [10101010, 123456789, 799999999, 990000000, 999999900];
    ais.bogus_mmsi = mmsi < 1e6 | ismember(mmsi, bogus_list);
end

% Flag if lon sign flip would make the path short
if has_tx && has_rx
    valid = ~isnan(ais.rxlat) & ~isnan(ais.rxlon) & ...
        ~isnan(ais.txlat) & ~isnan(ais.txlon);
    if any(valid)
        dist_flip = greatcircle(ais.rxlat(valid), ais.rxlon(valid), ...
            ais.txlat(valid), -ais.txlon(valid));
        flag_flip = false(nrows, 1);
        flag_flip(valid) = dist_flip < 1000;
        ais.likely_corrupt_lon = flag_flip;
    end
end
end

%% Helper to classify links
function ais = classify_links(ais, duct_strength_threshold)
nrows = height(ais);
var_names = ais.Properties.VariableNames;

cls = repmat("other", nrows, 1);

% Bad conditions override everything else
bad_fields = {'bad_receiver', 'bad_position', 'aircraft', ...
    'likely_corrupt_lon', 'known_bad_ship', 'bogus_mmsi'};
bad_mask = false(nrows, 1);
for i = 1:numel(bad_fields)
    if ismember(bad_fields{i}, var_names)
        bad_mask = bad_mask | logical(ais.(bad_fields{i}));
    end
end
cls(bad_mask) = "bad";

% Tropo if we have a duct strength value and not already bad
tropo_mask = ais.max_duct_strength > duct_strength_threshold; 
cls(~bad_mask & tropo_mask) = "tropo";

ais.link_classification = cls;
end

%% Calculate nmax implied by reflection angle
function ais = add_implied_nmax(ais, ais_freq_MHz)
ais.fmax_MHz = ais_freq_MHz * sind(abs(ais.elevation_midpoint_deg));
ais.nmax_elm3 = (ais.fmax_MHz * 1E6 / 9) .^ 2;
end

%% Helper to compute midpoint elevation angle for a skywave reflection
function ais = add_midpoint_elevation(ais, reflect_alt_km)

lat = ais.lat_mid;
lon = ais.lon_mid;
lat0 = ais.rxlat;
lon0 = ais.rxlon;

[~, elev, ~] = geodetic2aer(lat0, lon0, 0, lat, lon, reflect_alt_km, ...
    wgs84Ellipsoid('km'), 'degrees');
ais.elevation_midpoint_deg = elev;
end



%% Helper to pull midpoint columns with sensible fallbacks
function ais = get_midpoints(ais)
var_names = ais.Properties.VariableNames;


if all(ismember({'lat_mid', 'lon_mid'}, var_names))
    return
end

assert(all(ismember({'txlat', 'txlon', 'rxlat', 'rxlon'}, var_names)), ...
    'Missing tx/rx locations in AIS table.');
[mid_lat, mid_lon] = midpointLatLon(ais.txlat, ais.txlon, ais.rxlat, ais.rxlon);
ais.lat_mid = mid_lat;
ais.lon_mid = mid_lon;
end


%% Helper to gather receiver coords and times
function [rxlat, rxlon, timevec] = get_rx_and_time(tbl, default_time)
if isempty(tbl)
    rxlat = [];
    rxlon = [];
    timevec = [];
    return
end

var_names = tbl.Properties.VariableNames;
assert(all(ismember({'rxlat', 'rxlon'}, var_names)), 'Missing rxlat/rxlon');
rxlat = tbl.rxlat;
rxlon = tbl.rxlon;
timevec = get_time_vec(tbl);
if isempty(timevec)
    timevec = repmat(default_time, height(tbl), 1);
end
end

%% Helper to extract a time vector from a table
function tvec = get_time_vec(tbl)
tvec = [];
candidates = {'time', 'received_time_datenum_utc', 'received_time_utc', 'received_time'};
for i = 1:numel(candidates)
    if ismember(candidates{i}, tbl.Properties.VariableNames)
        tvec = tbl.(candidates{i});
        return
    end
end
end
