%% Specify filenames, times etc.
in_fname = 'data/ais/USN_JHUAPL_TER_%s_filtered_by_Area.nc';
times = datenum(2021, 7, 13):1/96:datenum(2021, 7, 16);
dist = 1000;

%% Load data
days = unique(round(times));
D = [];
D.time = [];
D.txlat= [];
D.txlon= [];
D.rxlat= [];
D.rxlon= [];

for d = 1:(length(days) - 1)
    in_fname_t = sprintf(in_fname, upper(datestr(days(d), 'ddmmmyyyy')));
    D.time = [D.time; ncread(in_fname_t, 'received_time_datenum_utc')];
    D.txlat = [D.txlat; ncread(in_fname_t, 'latitude_deg_n')];
    D.txlon = [D.txlon; ncread(in_fname_t, 'longitude_deg_e')];
    D.rxlat = [D.rxlat; ncread(in_fname_t, 'receiving_station_lat')];
    D.rxlon = [D.rxlon; ncread(in_fname_t, 'receiving_station_lon')];
end

%% calculate great circle distances
wgs84 = referenceEllipsoid('wgs84');
D.dist = distance(D.txlat, D.txlon, D.rxlat, D.rxlon, wgs84, 'degrees') / 1E3;


%% Histogram plot
% idx_list = {(D.dist < 1000); (D.dist > 1000)}; 
% clr = {'b', 'y'};
% labels = {'A (<1000 km)', 'B (>1000 km)'};
% clf
% hold on
% for i = 1:length(idx_list)
%     histogram(D.dist(idx_list{i}), 'binWidth', 20, 'FaceColor', clr{i});
%     set(gca,'YScale','log')
% end
% legend(labels)
% hold off
% set(gca, 'FontSize', 20)
% xlabel('Great Circle Distance (km)')
% ylabel('# 162-MHz AIS links')
% xlim(gca, [0, 2800])
% 
% grid on 
% 

%% remove nans and short links

goodind = D.dist > dist; % no >X links from NaN locations
fn = fields(D);
for i = 1:numel(fn)
    D2.(fn{i}) = D.(fn{i})(goodind);
end

%% time filter
links = zeros(1, length(times) - 1);
for t = 1:(length(times) - 1)

    tind = (D2.time >= times(t)) & (D2.time < times(t + 1));
    Dt = [];
    for i = 1:numel(fn)
        Dt.(fn{i}) = D2.(fn{i})(tind);
    end
    links(t) = length(Dt.dist);

end

%% Plot time series of link data
% clf
% ang = solarPosition(times(1:end-1), 44.4, -68.2,  0, 0, false);
% hold on
% 
% plot(times(1:end-1), links, '-xb', 'LineWidth', 3)
% for t = 1:size(ang, 1)
%     if ang(t, 1) > 90
%         plot([times(t), times(t)], [0 100000], '--k', 'LineWidth', 0.1)
%     end
% end
% 
% legend({sprintf('>%i-km AIS links', dist), 'Nighttime (Bar Harbor ME)'})
% hold off
% ylabel('Links identified in 15 mins')
% xlabel('UT (Date)')
% ylim([0, 800])
% datetick('x', 'HH (mmm. dd)', 'keeplimits')
% grid on
% grid minor
% set(gca, 'FontSize', 20)

%% plot links on globe
% idx = links == max(links(links < 600));
latS = -89:89;
lonS = -180:180;
sza = zeros(length(latS), length(lonS));
[latS_2d, lonS_2d] = meshgrid(latS, lonS);


for t = 1:length(times)
    %% Identify links
    time = times(t);
    disp(datestr(time))

    

    tind = (D2.time >= time) & (D2.time < time + 1/96);
    Dt = [];
    for i = 1:numel(fn)
        Dt.(fn{i}) = D2.(fn{i})(tind);
    end

    latI_f = [];
    lonI_f = [];
    for s = 1:length(Dt.time)

        lat = [Dt.txlat(s), Dt.rxlat(s)];
        lon = [Dt.txlon(s), Dt.rxlon(s)];
        [latI,lonI] = interpm(lat,lon,1,'gc');
        latI_f = [latI_f; latI; NaN];
        lonI_f = [lonI_f; lonI; NaN];
    end

    %% Plot links
    uif = uifigure;
    g = geoglobe(uif);
    [camlat,camlon,camh] = campos(g); 
    geoplot3(g, latI_f, lonI_f, ones(size(latI_f)) * 2000,'y', 'LineWidth', 2)

    %% Plot the solar terminator
    for i = 1:length(latS)
        for j = 1:length(lonS)
            ang = solarPosition(time, latS(i), lonS(j), 0, 0, false);
            sza(i, j) = ang(1);
        end
    end
    geoplot3(g, latS_2d(sza'>90), lonS_2d(sza'>90), ones([sum(sza(:)>90), 1]) * 1000, 'ko')

    %% Set the camera
    campos(g,30,-70,camh);

    %% Save to PDF
    exportapp(uif, filename('plots/{yyyymmmdd_HHMM}.pdf', time))

    close(findall(0, 'type', 'figure'));
    
end















