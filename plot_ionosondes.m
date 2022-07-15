%%
times = datenum(2021, 7, 13:15);
fname_fmt = 'data/ionosonde/{NAME}_{yyyymmdd}.nc';
sites = {'MILLSTONEHILL', 'WALLOPSIS', 'RAMEY'};
vars = {'foEs', 'Time'};


%% load
D = [];
for s = 1:length(sites)
    for v = vars
        D{s}.(v{1}) = [];
    end
    for t = 1:length(times)
        fname = filename(fname_fmt, times(t), sites{s});
        for v = vars

            D{s}.(v{1}) = [D{s}.(v{1}); ncread(fname, v{1})];

        end
        D{s}.lat = ncreadatt(fname, '/', 'lat');
        D{s}.lon= ncreadatt(fname, '/', 'lon');
    end
end

%% SZA




%% plot ionosonde foEs data
clf
for s = 1:length(sites)
    subplot(length(sites), 1, s)
    ang = solarPosition(D{s}.Time, D{s}.lat, D{s}.lon,  0, 0, false);

    hold on
    plot(D{s}.Time, D{s}.foEs, '-xb', 'LineWidth', 2)
    for t = 1:size(ang, 1)
        if ang(t, 1) > 90
            plot([D{s}.Time(t), D{s}.Time(t)], [0 12], '--k', 'LineWidth', 0.1)
        end
    end
    ylim([0, 12])
    ylabel(sprintf('%s\nFoEs (MHz) ', sites{s}))
    if s < length(sites)
        set(gca, 'XTickLabels', '')
        
    else
        datetick('x', 'HH (mmm. dd)', 'keeplimits')
        xlabel('UT (Date)')

    end
        legend({"FoEs", "Nighttime"})
    grid on
    grid minor
    set(gca, 'FontSize', 20)
    hold off
end


%% plot locations on globe

latI_f = [];
lonI_f = [];
for s = 1:length(sites)

    latI_f = [latI_f, D{s}.lat];
    lonI_f = [lonI_f,  D{s}.lon];
end
uif = uifigure;
g = geoglobe(uif);

%%
geoplot(latI_f, lonI_f, ones(size(latI_f)) * 2000,'ro', MarkerSize=20)

