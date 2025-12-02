% Script to export hourly ERA5 maximum duct strength (M-excess) grids
% into individual NetCDF files. Follows the usage pattern from
% example_load_and_plot_era5_mexcess.m.

%% Set inputs
start_time = datenum([2024 5 8 0 0 0]);
end_time = datenum([2024 5 13 0 0 0]);
lat_box = [10 80]; % deg N
lon_box = [-180 -50]; % deg E
outputFolder = '~/data/sporadice/era5/';
overwriteExistingFiles = true; % default: no-clobber, skip existing files


latlonType = 'abs';
modelGrid = 'era5_f320';
nPoints = 1;
varNums = [114 116];
datatype = 0;
precisionFlag = 2;
xyThin = 1;
nzRead = inf;
readAllData = 0;
rootFileFolder = '/Volumes/era5';


hourly_times = start_time:(1/24):(end_time - 1/24);
nHours = numel(hourly_times);

%%  Output configuration
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

for idx = 1:nHours
    this_time_datenum_utc = hourly_times(idx);
    timeLabel = datestr(this_time_datenum_utc, 'yyyy-mm-dd HH:MM');
    fprintf('(%d/%d) Processing %s UTC\n', idx, nHours, timeLabel);

    outputFilename = fullfile(outputFolder, ...
        sprintf('max_duct_strength_%s.nc', datestr(this_time_datenum_utc, 'yyyymmdd_HHMM')));

    if exist(outputFilename, 'file')
        if overwriteExistingFiles
            delete(outputFilename);
        else
            fprintf('    %s already exists, skipping.\n', outputFilename);
            continue;
        end
    end

    tempFilename = [outputFilename '.tmp'];
    if exist(tempFilename, 'file')
        delete(tempFilename);
    end

    startTimer = tic;
    try
        % Pull ERA5 data exactly as in the plotting example
        D = get_nwp_data_newh5(lat_box, lon_box, latlonType, modelGrid, nPoints, ...
            varNums, datatype, this_time_datenum_utc, precisionFlag, xyThin, ...
            nzRead, readAllData, rootFileFolder);

        if min(D.lonvec) > 180
            D.lonvec = D.lonvec - 360;
        end

        % Replace non-duct points (height == 0) with NaNs before squeezing
        D.data.m_excess_ref15m_ha(D.data.ducthgt_ha == 0) = nan;
        max_duct_strength = squeeze(D.data.m_excess_ref15m_ha);
        duct_height_km = single(D.data.ducthgt_ha);
        duct_height_km(duct_height_km == 0) = nan;
        duct_height_km = duct_height_km ./ single(1000); % convert meters -> kilometers

        % Create the NetCDF file and write coordinates + data via temp file
        nccreate(tempFilename, 'lat', 'Dimensions', {'lat', length(D.latvec)}, ...
            'Datatype', 'double', 'Format', 'netcdf4');
        ncwrite(tempFilename, 'lat', D.latvec);
        ncwriteatt(tempFilename, 'lat', 'units', 'degrees_north');

        nccreate(tempFilename, 'lon', 'Dimensions', {'lon', length(D.lonvec)}, ...
            'Datatype', 'double');
        ncwrite(tempFilename, 'lon', D.lonvec);
        ncwriteatt(tempFilename, 'lon', 'units', 'degrees_east');

        nccreate(tempFilename, 'time', 'Dimensions', {'time', 1}, ...
            'Datatype', 'double');
        hours_since_1900 = (this_time_datenum_utc - datenum(1900, 1, 1)) * 24;
        ncwrite(tempFilename, 'time', hours_since_1900);
        ncwriteatt(tempFilename, 'time', 'units', 'hours since 1900-01-01 00:00:00');
        ncwriteatt(tempFilename, 'time', 'standard_name', 'time');

        nccreate(tempFilename, 'max_duct_strength', ...
            'Dimensions', {'lat', length(D.latvec), 'lon', length(D.lonvec)}, ...
            'Datatype', 'single', 'FillValue', single(1e20));
        ncwrite(tempFilename, 'max_duct_strength', single(max_duct_strength));
        ncwriteatt(tempFilename, 'max_duct_strength', 'units', 'dB');
        ncwriteatt(tempFilename, 'max_duct_strength', 'long_name', ...
            'Maximum duct strength (M-excess) at 15 m reference height');
        
        nccreate(tempFilename, 'duct_height', ...
            'Dimensions', {'lat', length(D.latvec), 'lon', length(D.lonvec)}, ...
            'Datatype', 'single', 'FillValue', single(1e20));
        ncwrite(tempFilename, 'duct_height', duct_height_km);
        ncwriteatt(tempFilename, 'duct_height', 'units', 'km');
        ncwriteatt(tempFilename, 'duct_height', 'long_name', 'Duct height');

        % Global attributes for provenance
        ncwriteatt(tempFilename, '/', 'title', ...
            'ERA5 maximum duct strength (M-excess)');
        ncwriteatt(tempFilename, '/', 'source', 'get_nwp_data_newh5 via save_hourly_era5_mexcess_netcdf.m');
        ncwriteatt(tempFilename, '/', 'history', ...
            ['Created on ' datestr(now, 'yyyy-mm-dd HH:MM') ' UTC']);
        ncwriteatt(tempFilename, '/', 'time_coverage_start', ...
            datestr(this_time_datenum_utc, 'yyyy-mm-ddTHH:MM:SSZ'));
        ncwriteatt(tempFilename, '/', 'geospatial_lat_min', lat_box(1));
        ncwriteatt(tempFilename, '/', 'geospatial_lat_max', lat_box(2));
        ncwriteatt(tempFilename, '/', 'geospatial_lon_min', lon_box(1));
        ncwriteatt(tempFilename, '/', 'geospatial_lon_max', lon_box(2));

        movefile(tempFilename, outputFilename);
        fprintf('    Completed %s in %.1f seconds.\n', timeLabel, toc(startTimer));
    catch ME
        fprintf('    Error processing %s after %.1f seconds: %s\n', ...
            timeLabel, toc(startTimer), ME.message);
        if exist(tempFilename, 'file')
            delete(tempFilename);
        end
        continue;
    end
end

fprintf('Done! NetCDF files located in %s\n', outputFolder);
