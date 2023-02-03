function plot_STATS_5minBefore5to10minafter_anova(cells1,cells2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats comparing Firing Rate over Time 5 min before injection and 5 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));

[cntrlDiff1,ketDiff1] = calc_FRDiff_5minBefore5to10minAfter(cells1);
[cntrlDiff2,ketDiff2] = calc_FRDiff_5minBefore5to10minAfter(cells2);


%%
y = vertcat(cntrlDiff1, ketDiff1, cntrlDiff2, ketDiff2);

x = vertcat(...
    repmat({'Control_excitatory'}, size(cntrlDiff1)),...
    repmat({'Ketamine_excitatory'}, size(ketDiff1)),...
    repmat({'Control_interneurons'}, size(cntrlDiff2)),...
    repmat({'Ketamine_interneurons'}, size(ketDiff2))...
    );
%%

close all; clear g;
figure('Position',[100 100 1200 600]);
g = gramm('y',y,'x',x,'color',x);
% g.stat_ellipse
g.stat_boxplot('notch',1);
% g.stat_violin('normalization','area','width',0.2,'extra_y',0)
g.draw();
%%
close all; clear g;
[p,t,stats] = anova1(y,x);
[results,~,~,gnames] = multcompare(stats,'CriticalValueType','bonferroni');
end

