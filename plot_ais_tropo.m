%% plot_ais_tropo.m
% Load AIS data

%%
times = datenum(2021, 7, 13):1/24:datenum(2021, 7, 15, 23, 0, 0);
lat = 44.38;
lon = -68.20;
fname = 'data/ais_tropo/eastcoast_{YYYYmmdd}_{HH}Z.nc';
tropo_ductstrength = zeros(size(times));

for t = 1:length(times)
    ds = ncread(filename(fname, times(t)), 'era5_duct_strength_m_units');
    tropo_ductstrength(t) = median(ds(ds>0));

end


%% plot
ang = solarPosition(times, lat, lon,  0, 0, false);
hold on
for t = 1:size(ang, 1)
    if isnan(ang(t, 2))
        plot([times(t), times(t)], [0 6], '--k', 'LineWidth', 0.1)
    end
end
plot(times, tropo_ductstrength, '-x', 'LineWidth', 3);
hold off
ylim([0, 6])
datetick('x', "mm/dd")
xlabel('Time (UT)')
ylabel('Median Tropo Duct Strength (M units)')
grid on
grid minor
set(gca, 'FontSize', 20)


%% 2d plots
time = datenum(2021, 7, 14, 19, 0, 0);
%time = datenum(2021, 7, 13, 14, 0, 0);
climit = [0, 20];
ds = ncread(filename(fname, time), 'era5_duct_strength_m_units');
lats = ncread(filename(fname, time), 'era5_latitude_deg_n');
lons = ncread(filename(fname, time), 'era5_longitude_deg_e');
ds(ds == 0) = nan;

m_proj('Miller Cylindrical', 'lat', [25, 50], 'long', [-85, -50]);

[~, hC] = m_contourf(lons, lats, squeeze(ds), 100);
set(hC, 'LineStyle', 'none')
m_coast('color', 'k', 'LineWidth', 1)
m_grid
caxis(climit)


figure
caxis(climit)
h = colorbar;
ylabel(h, 'Tropo Duct Strength (M Units)')
sfh3.Visible = 'off';
set(gca, 'FontSize', 26, 'Fontweight', 'bold')



























