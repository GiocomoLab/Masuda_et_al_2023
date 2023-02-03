function plot_STATS_comparing_wt_hcn1ko_mk801_before_and_after_dchPeriod(wt_ket_Cells,hcn1ko_ket_Cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats comparing Firing Rate over during decoherence period and
% during equivalent length of time before dch period for wt, hcn1ko, and
% mk801
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));

dchAndBefore = find_dch_and_equivalent_beforetime(wt_ket_Cells_onlyStable);
wt_ket_change_in_dchFR =  dchAndBefore.dch_avg_cellFR - dchAndBefore.before_dch_avg_cellFR;


dchAndBefore = find_dch_and_equivalent_beforetime(hcn1ko_ket_Cells_onlyStable);
hcn1ko_ket_change_in_dchFR =  dchAndBefore.dch_avg_cellFR - dchAndBefore.before_dch_avg_cellFR;

% dchAndBefore = find_dch_and_equivalent_beforetime(wt_mk801_Cells);
% wt_mk801_change_in_dchFR =  dchAndBefore.dch_avg_cellFR - dchAndBefore.before_dch_avg_cellFR;
% 
% dchAndBefore = find_dch_and_equivalent_beforetime(hcn1ko_mk801_Cells);
% hcn1ko_mk801_change_in_dchFR =  dchAndBefore.dch_avg_cellFR - dchAndBefore.before_dch_avg_cellFR;

color = horzcat(...
    repmat({'wt_ket'}, size(wt_ket_change_in_dchFR)),...
    repmat({'hcn1ko_ket'}, size(hcn1ko_ket_change_in_dchFR))...
%     repmat({'wt_mk801'}, size(wt_mk801_change_in_dchFR)),...
%     repmat({'hcn1ko_mk801'}, size(hcn1ko_mk801_change_in_dchFR))...
    );
% y = horzcat(wt_ket_change_in_dchFR,hcn1ko_ket_change_in_dchFR,wt_mk801_change_in_dchFR,hcn1ko_mk801_change_in_dchFR);
y = horzcat(wt_ket_change_in_dchFR,hcn1ko_ket_change_in_dchFR);
%%
close all;
clear g;
g(1,1) = gramm('x',color,'y',y,'color',color);
g(1,1).stat_boxplot();
g.draw()

[p,h] = ranksum(wt_ket_change_in_dchFR,hcn1ko_ket_change_in_dchFR)
% [p,h] = ranksum(wt_mk801_change_in_dchFR,hcn1ko_mk801_change_in_dchFR)
% [p,h] = ranksum(wt_ket_change_in_dchFR,wt_mk801_change_in_dchFR)
% [p,h] = ranksum(hcn1ko_ket_change_in_dchFR,hcn1ko_mk801_change_in_dchFR)