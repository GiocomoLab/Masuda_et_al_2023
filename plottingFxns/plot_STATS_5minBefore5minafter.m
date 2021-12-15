function plot_STATS_5minBefore5minafter(cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats comparing Firing Rate over Time 5 min before injection and 5 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
sampleRate = 50; %hz
secInMin = 60; 
scaling  = sampleRate * secInMin; 
% dsfactor = 50;
%%
cntrl = cells.timeFRcircaControlInjx;
ket = cells.timeFRcircaKetamineInjx;

cntrlBefore = cntrl(:,1:scaling*5);
ketBefore = ket(:,1:scaling*5);

cntrlAfter = cntrl(:,scaling * 5+1:scaling * 10);
ketAfter = ket(:,scaling * 5+1:scaling * 10);

% cntrlDiff = mean(cntrlAfter,2)-mean(cntrlBefore,2);
% ketDiff = mean(ketAfter,2)-mean(ketBefore,2);

cntrlDiff = mean(cntrlAfter-cntrlBefore,2);
ketDiff = mean(ketAfter-ketBefore,2);
% 
% cntrlDiff = mean((cntrlAfter-cntrlBefore)/nanmean(cntrlBefore),2);
% ketDiff = mean((ketAfter-ketBefore)/nanmean(cntrlBefore),2);

% absCntrlDiff = abs(cntrlDiff);
% absKetDiff = abs(cntrlDiff);

[h,p] = ttest2(cntrlDiff, ketDiff);
fprintf(num2str(p))
p = ranksum(cntrlDiff,ketDiff);
fprintf('\n')
fprintf(num2str(p))

ket_sem = std( ketDiff ) / sqrt( numel(ketDiff) );
cntrl_sem = std( cntrlDiff ) / sqrt( numel(cntrlDiff) );
%%
close all; clear g;
figure(); hold on;
y = vertcat(cntrlDiff,ketDiff);

x = vertcat(...
    repmat({'Control'}, size(cntrlDiff)),...
    repmat({'Ketamine'}, size(ketDiff))...
    );

boxplot(y,x,'ColorGroup',x,'colors','k','symbol','','PlotStyle','traditional','Widths',0.5)
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
ylabel('Hz')
set(findobj(gca,'type','line'),'linew',2)
title(sprintf('P-value: %0.9f',p))
end

