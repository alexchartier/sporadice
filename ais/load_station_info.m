%% load station info lat/lon
function station_info = load_station_info(station_fn)
station_info = [];

if isempty(station_fn) || ~isfile(station_fn)
    fprintf('    Station info file missing: %s\n', station_fn);
    return
end

raw = load(station_fn);
fn = fieldnames(raw);

if ismember('lat', fn) && ismember('lon', fn)
    station_info.lat = raw.lat;
    station_info.lon = raw.lon;
    return
end

if numel(fn) == 1 && isstruct(raw.(fn{1})) && ...
        isfield(raw.(fn{1}), 'lat') && isfield(raw.(fn{1}), 'lon')
    station_info.lat = raw.(fn{1}).lat;
    station_info.lon = raw.(fn{1}).lon;
    return
end

error('Station info in %s missing lat/lon fields.', station_fn);
