%% Plot links
for t = 1:length(times)
    x = load(filename(out_fn_fmt, times(t)));
    Tbl_far_good = x.Tbl_far_good;
    tidx = abs(Tbl_far_good.time - times(t)) < 15/60/24;
    C = load('coastlines');
    nans = nan(size(Tbl_far_good.rxlat(tidx)));
    latarr = [Tbl_far_good.txlat(tidx), Tbl_far_good.rxlat(tidx), nans]';
    lonarr = [Tbl_far_good.txlon(tidx), Tbl_far_good.rxlon(tidx), nans]';
    hold on
    plot(cell2mat(station_info(:, 3)), cell2mat(station_info(:, 2)),...
        '.k', 'MarkerSize',10)
    plot(lonarr(:), latarr(:), '-r')
    plot(C.coastlon, C.coastlat,'k')


    title(datestr(times(t), 'yyyy/mm/dd HH:MM'))
    hold off
    xlim([-180, -30])
    ylim([0, 60])
    
    pause(0.6)
    clf
end