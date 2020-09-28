function [H,P] = plot_dch_autocorrelationScore_cells1_VS_cells2(cells1,cells2,name1,name2)


addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
%%

cells1_autocorrelationScore = calc_dch_autocorrelationScore(cells1);
cells2_autocorrelationScore = calc_dch_autocorrelationScore(cells2);

[H,P] = ttest2(cells1_autocorrelationScore,cells2_autocorrelationScore);

%%
close all; clear g;
y = vertcat(cells1_autocorrelationScore,cells2_autocorrelationScore);

x = vertcat(...
    repmat({name1}, size(cells1_autocorrelationScore)),...
    repmat({name2}, size(cells2_autocorrelationScore))...
    );

g = gramm('x',y,'color',x);
g.stat_bin('geom','stairs','fill','transparent');
g.set_title(sprintf('P-value: %0.3f',P));
g.draw;
end