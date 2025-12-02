%% Set inputs
% in_fn_fmt = '~/data/ais/{yyyymm}/{dd}/uscg_nais_{yyyymmdd_HH}.mat';
in_fn_fmt = '/Volumes/processed/{yyyymm}/{dd}/uscg_nais_{yyyymmdd_HH}.mat';
out_fn_fmt = '~/data/sporadice/ais_proc_l1/{yyyymmdd_HHMM}.mat';
station_fn = '~/data/sporadice/ais_station_info/station_info.mat';
times = datenum(2022, 5, 1):1/24:datenum(2025, 5, 1)-1E-9;
far = 1E3; % How many kilometres is 'far' away

%% Load
for t = 1:length(times)
    disp(datestr(times(t)))
    D = loadstruct(filename(in_fn_fmt, times(t)));
    ais = D.ais;
    clear Dan
    station_info = ais.info.nais_station_info;


    %% Calculate great circle distances
    Tbl = table;
    Tbl.time = ais.received_time_datenum_utc;
    Tbl.mmsi = ais.mmsi;
    Tbl.txlat = ais.latitude_deg_n;
    Tbl.txlon = ais.longitude_deg_e;
    Tbl.rxlat = cell2mat(station_info(ais.receiving_station_index,2));
    Tbl.rxlon = cell2mat(station_info(ais.receiving_station_index,3));

    Tbl.gcdist = greatcircle(Tbl.rxlat, Tbl.rxlon, ...
        Tbl.txlat, Tbl.txlon);

    %% Check the far-away ships are near where we'd expect them to be
    far_idx = Tbl.gcdist > far;
    Tbl_far = Tbl(far_idx, :);
    % badid = Tbl_far.txlat == 0 | Tbl_far.txlon == 0;
    % Tbl_far = Tbl_far(~badid, :);

    far_ships = unique(Tbl_far.mmsi);

    % Calculate the median location of each far ship from the full table
    Tbl_far.txlat_med = nan(size(Tbl_far.txlat));
    Tbl_far.txlon_med = nan(size(Tbl_far.txlon));

    for i = 1:length(far_ships)
        ship_idx = Tbl.mmsi == far_ships(i);
        far_ship_idx = Tbl_far.mmsi == far_ships(i);

        % Require at least 4 measurements of the ship
        if sum(ship_idx) < 3
            continue
        end

        % Require more than one location in the full table (anti-spoofing)
        if isscalar(unique(Tbl(ship_idx, 'txlat').txlat)) || ...
                isscalar(unique(Tbl(ship_idx, 'txlon').txlon))
            continue
        end
        Tbl_far.txlat_med(far_ship_idx) = median(Tbl(ship_idx, 'txlat').txlat);
        Tbl_far.txlon_med(far_ship_idx) = median(Tbl(ship_idx, 'txlon').txlon);

    end

    % Check great-circle distances of 'far' receptions are not too far from the
    % median locations
    gcdist_check = greatcircle(Tbl_far.txlat, Tbl_far.txlon, ...
        Tbl_far.txlat_med, Tbl_far.txlon_med);
    Tbl_far_good = Tbl_far(gcdist_check < far, :);

    %% Save
    save(filename(out_fn_fmt, times(t)), "Tbl_far_good")

end

%% Save out the station information
out.lat = round(cell2mat(station_info(:, 2)) * 1E3) / 1E3;
out.lon = round(cell2mat(station_info(:, 3)) * 1E3) / 1E3;

finind = ~isnan(out.lat);
out.lat = out.lat(finind);
out.lon = out.lon(finind);
savestruct(station_fn, out)




