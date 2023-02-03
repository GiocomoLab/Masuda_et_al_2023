function [H,P] = plot_dch_autocorrelationScore_cells1_VS_cells2(cells1,cells2,name1,name2)


addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
%%

cells1_autocorrelationScore = calc_dch_autocorrelationScore_old(cells1);
cells2_autocorrelationScore = calc_dch_autocorrelationScore_old(cells2);
% [cells1_autocorrelationScore,p1_cntrl, p1_dch] = calc_dch_autocorrelationScore(cells1);
% [cells2_autocorrelationScore,p2_cntrl, p2_dch] = calc_dch_autocorrelationScore(cells2);
[H,P] = ttest2(cells1_autocorrelationScore,cells2_autocorrelationScore);

%%
figure();
clear g;
y = vertcat(cells1_autocorrelationScore,cells2_autocorrelationScore);

x = vertcat(...
    repmat({name1}, size(cells1_autocorrelationScore)),...
    repmat({name2}, size(cells2_autocorrelationScore))...
    );

g = gramm('x',y,'color',x);
g.stat_bin('geom','bar','fill','transparent');
g.set_title(sprintf('P-value: %0.3f',P));
g.draw;
%%
% figure();
% clear g;
% y = vertcat(p1_cntrl,p2_cntrl);
% 
% x = vertcat(...
%     repmat({name1}, size(p1_cntrl)),...
%     repmat({name2}, size(p2_cntrl))...
%     );
% [H,P] = ttest2(p1_cntrl,p2_cntrl);
% g = gramm('x',x,'y',y);
% g.stat_violin('normalization','width','dodge',0,'fill','edge');
% g.stat_boxplot('width',0.15);
% g.set_title(sprintf('Cntrl P-value: %0.3f',P));
% g.draw;
% %%
% figure();
% clear g;
% y = vertcat(p1_dch,p2_dch);
% 
% x = vertcat(...
%     repmat({name1}, size(p1_dch)),...
%     repmat({name2}, size(p2_dch))...
%     );
% [H,P] = ttest2(p1_dch,p2_dch);
% g = gramm('x',x,'y',y);
% g.stat_violin('normalization','width','dodge',0,'fill','edge');
% g.stat_boxplot('width',0.15);
% g.set_title(sprintf('Dch P-value: %0.3f',P));
% g.draw;

%%
figure();
clear g;
y = cells1_autocorrelationScore;

x = repmat({name1}, size(cells1_autocorrelationScore));

g = gramm('x',y,'color',x);
g.stat_bin('geom','bar','fill','transparent');
g.set_color_options('chroma',0);
g.set_title('Change in Peak Prominence (+ is more decohered)');
g.draw;
%%
end