close all
clear all

% example script to pull ERA5 M-excess data and overlay on a map
% written by T.R. Hanley, last updated 4-Nov-2025

% Lat/Lon Box that covers all of NAIS areas with exception of Guam 
lat_box = [10 80]; % deg N
lon_box = [-180 -50]; % deg E %Note that this doesn't work crossing 180E/180W


this_time_datenum_utc = datenum([2024 5 10 0 0 0]); % yyyy mm dd HH MM SS

latlonType = 'abs';
modelGrid = 'era5_f320';
nPoints = 1;
varNums = [114 116];
datatype = 0;
precisionFlag = 2;
xyThin = 1;
nzRead = inf;
readAllData = 0;
rootFileFolder = '\\isilonshares\AMDS\A2\A2A\Hanley\epc_storage\metoc\nwp\reanalysis\era5'; % replace with appropriate path

D = get_nwp_data_newh5( lat_box, lon_box, latlonType, modelGrid, nPoints,...
    varNums, datatype, this_time_datenum_utc, ...
    precisionFlag, xyThin, nzRead, ...
    readAllData, rootFileFolder );

if min(D.lonvec)>180
    D.lonvec = D.lonvec - 360;
end
D.data.m_excess_ref15m_ha(D.data.ducthgt_ha==0) = nan; % if no duct is present (based on a height of 0), set to NaN

cols = wide_range_colormap(24);

max_duct_strength = squeeze(D.data.m_excess_ref15m_ha);
max_duct_strength(max_duct_strength==0) = -3; % set for colorbar scaling so values close to zero gets assigned color below zero

hf = figure;
imagesc(D.lonvec,D.latvec,max_duct_strength);
set(gca,'ydir','normal');

colormap([1 1 1; cols]);
caxis([-30 95]);
ylim(lat_box);
add_coast;
hold on;
grid on;
ylabel(['Latitude (' char(176) 'N)'],'FontSize',14);
xlabel(['Longitude (' char(176) 'E)'],'FontSize',14);

ha = gca;
set(ha,'FontSize',14);
hc = colorbar;
set(hc,'FontSize',14,'ytick',[-25:5:95]);
hcl = get(hc,'ylabel');
set(hcl,'string','M-Excess ( \Delta M)','FontSize',14)
set(ha,'xlim',lon_box);
title({['M-Excess: ' datestr(this_time_datenum_utc,'yyyy-mmm-dd HH:MM') ' UTC']; 'From ERA5 Weather Model'},'FontSize',16);
