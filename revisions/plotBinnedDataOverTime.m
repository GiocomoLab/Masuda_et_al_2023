function plotBinnedDataOverTime(x, ctrl, ket, varargin)

figure; clf; hold on;
if nargin>3
    base = varargin{2};
    %base(base>20) = NaN;
    data_SEM = nanstd(base,1)./sqrt(size(base,1));
    shadedErrorBar(x,nanmean(base,1),data_SEM,'lineprops','-k');
end
%ctrl(ctrl>20) = NaN;
data_SEM = nanstd(ctrl,1)./sqrt(size(ctrl,1));
shadedErrorBar(x,nanmean(ctrl,1),data_SEM,'lineprops','-m');
%ket(ket>20) = NaN;
data_SEM = nanstd(ket,1)./sqrt(size(ket,1));
shadedErrorBar(x,nanmean(ket,1),data_SEM,'lineprops','-g');
vline(0,'k');

% smoothSigma = 5;
% smoothWindow = floor(smoothSigma*5/2)*2+1;
% gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
% rescale = 60;
% x = 1:30;
% for a = 1:size(ket,1)
%     if nargin>3
%         base_smoothed(a,:) = nanmean(reshape(base(a,1:rescale*fix(length(base)/rescale)),rescale,[]));%conv(base(a,:),gauss_filter,'same');
%     end
%     ctrl_smoothed(a,:) = nanmean(reshape(ctrl(a,1:rescale*fix(length(ctrl)/rescale)),rescale,[]));%conv(ctrl(a,:),gauss_filter,'same');
%     ket_smoothed(a,:) = nanmean(reshape(ket(a,1:rescale*fix(length(ket)/rescale)),rescale,[]));%conv(ket(a,:),gauss_filter,'same');
% end
% 
% figure; clf; hold on;
% if nargin>3
%     data_SEM = nanstd(base_smoothed,1)./sqrt(size(base_smoothed,1));
%     shadedErrorBar(x,nanmean(base_smoothed,1),data_SEM,'lineprops','-k');
% end
% data_SEM = nanstd(ctrl_smoothed,1)./sqrt(size(ctrl_smoothed,1));
% shadedErrorBar(x,nanmean(ctrl_smoothed,1),data_SEM,'lineprops','-m');
% data_SEM = nanstd(ket_smoothed,1)./sqrt(size(ket_smoothed,1));
% shadedErrorBar(x,nanmean(ket_smoothed,1),data_SEM,'lineprops','-g');
% vline(0,'k');


