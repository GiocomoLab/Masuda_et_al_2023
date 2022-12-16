function plot_peakFR(cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats comparing Firing Rate over Time 5 min before injection and 5 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
sampleRate = 50; %hz
secInMin = 60; 
scaling  = sampleRate * secInMin; 

% dsfactor = 50;
%
cntrl = cells.timeFRcircaControlInjx;
ket = cells.timeFRcircaKetamineInjx;

% cntrlAfter = cntrl(:,scaling * 5+1:scaling * 10);
ketAfter = ket(:,scaling * 5+1:scaling * 10);
numCells = size(smooth_ketAfter,1);
smooth_ketAfter = smoothdata(ketAfter,2,"movmean",scaling);
plot(smooth_ketAfter(1,:)')
peakTimes = nan(numCells,1);
for i = 1:numCells
    [peaks,locs,~,~] = findpeaks(smooth_ketAfter(i,:),'MinPeakProminence',2);
    [~,I] = max(peaks);
    if isempty(I)
        continue
    end
    peakTimes(i) = locs(I);
end
meanPeakTime = nanmean(peakTimes./scaling);
stdPeakTime = nanstd(peakTimes./scaling);
stderror = stdPeakTime / sqrt( numCells );
%%
% close all; clear g;
% figure(); hold on;
% y = vertcat(cntrlDiff,ketDiff);
% 
% x = vertcat(...
%     repmat({'Control'}, size(cntrlDiff)),...
%     repmat({'Ketamine'}, size(ketDiff))...
%     );
% 
% boxplot(y,x,'ColorGroup',x,'colors','k','symbol','','PlotStyle','traditional','Widths',0.5)
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% ylabel('Hz')
% set(findobj(gca,'type','line'),'linew',2)
% title(sprintf('P-value: %0.9f',p))
end

