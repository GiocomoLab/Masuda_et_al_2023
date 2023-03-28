function plotBinnedDataOverTime(x, ctrl, ket, fig_title)

figure; clf; hold on;
data_SEM = nanstd(ctrl,1)./sqrt(size(ctrl,1));
shadedErrorBar(x,nanmean(ctrl,1),data_SEM,'lineprops','-m');
data_SEM = nanstd(ket,1)./sqrt(size(ket,1));
shadedErrorBar(x,nanmean(ket,1),data_SEM,'lineprops','-g');
vline(0,'k');
title(fig_title)

end