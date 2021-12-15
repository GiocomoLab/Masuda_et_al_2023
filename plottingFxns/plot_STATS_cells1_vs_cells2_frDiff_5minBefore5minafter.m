function plot_STATS_cells1_vs_cells2_frDiff_5minBefore5minafter(cells1,label1,cells2,label2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats comparing Firing Rate over Time 5 min before injection and 5
% min after injection between two subpopulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));

[cntrlDiff1,ketDiff1] = find_frDiff_5minBefore5minafter(cells1);
[cntrlDiff2,ketDiff2] = find_frDiff_5minBefore5minafter(cells2);

[h,p] = ttest2(ketDiff1, ketDiff2);
fprintf(num2str(p))
p = ranksum(ketDiff1,ketDiff2);
fprintf('\n')
fprintf(num2str(p))



ket_sem1 = std( ketDiff1 ) / sqrt( numel(ketDiff1) );
ket_sem2 = std( ketDiff2 ) / sqrt( numel(ketDiff2) );
cntrl_sem1 = std( cntrlDiff1 ) / sqrt( numel(cntrlDiff1) );
cntrl_sem2 = std( cntrlDiff2 ) / sqrt( numel(cntrlDiff2) );
%%
close all; clear g;
figure(); hold on;
y = vertcat(cntrlDiff1,cntrlDiff2,ketDiff1,ketDiff2);

x = vertcat(...
    repmat({strcat(label1,'_Cntrl')}, size(cntrlDiff1)),...
    repmat({strcat(label2,'_Cntrl')}, size(cntrlDiff2)),...
    repmat({strcat(label1,'_Ket')}, size(ketDiff1)),...
    repmat({strcat(label2,'_Ket')}, size(ketDiff2))...
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
title(sprintf('P-value: %0.2f',p))
%% Anova Stats
drug = vertcat(...
    repmat({strcat('Cntrl')}, size(cntrlDiff1)),...
    repmat({strcat('Cntrl')}, size(cntrlDiff2)),...
    repmat({strcat('Ket')}, size(ketDiff1)),...
    repmat({strcat('Ket')}, size(ketDiff2))...
    );


celltype = vertcat(...
    repmat({label1}, size(cntrlDiff1)),...
    repmat({label2}, size(cntrlDiff2)),...
    repmat({label1}, size(ketDiff1)),...
    repmat({label2}, size(ketDiff2))...
    );
[~,~,stats] = anovan(y,{drug,celltype},'model','interaction',...
    'varnames',{'drug','celltype'});

[results,means] = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
end

