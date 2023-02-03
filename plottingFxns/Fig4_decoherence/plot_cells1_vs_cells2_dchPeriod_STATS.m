function plot_cells1_vs_cells2_dchPeriod_STATS(cells1,label1,cells2,label2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats comparing Firing Rate over during decoherence period and during equivalent length of time before dch period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
dchAndBefore1 = find_dch_and_equivalent_beforetime(cells1);
dchAndBefore2 = find_dch_and_equivalent_beforetime(cells2);
%%
% abs_fr_diff1 = abs(dchAndBefore1.dch_avg_cellFR - dchAndBefore1.before_dch_avg_cellFR);
% abs_fr_diff2 = abs(dchAndBefore2.dch_avg_cellFR - dchAndBefore2.before_dch_avg_cellFR);

fr_diff1 = (dchAndBefore1.dch_avg_cellFR - dchAndBefore1.before_dch_avg_cellFR)./nanmean(dchAndBefore1.before_dch_avg_cellFR);
fr_diff2 = (dchAndBefore2.dch_avg_cellFR - dchAndBefore2.before_dch_avg_cellFR)./nanmean(dchAndBefore2.before_dch_avg_cellFR);


nanIndx1 = isnan(fr_diff1);
fr_diff1 = fr_diff1(~nanIndx1);
% metadata1 = cells1.metadata(:,4);
% metadata1 = metadata1(~nanIndx1);
metadata1 = repmat({label1},size(fr_diff1));

nanIndx2 = isnan(fr_diff2);
fr_diff2 = fr_diff2(~nanIndx2);
% metadata2 = cells2.metadata(:,4);
% metadata2 = metadata2(~nanIndx2);
metadata2 = repmat({label2},size(fr_diff2));

[p,h] = ranksum(fr_diff1, fr_diff2);
fprintf('The hypothesis test is %i; P-value is %.02d\n',h,p);


fr_diff = horzcat(fr_diff1,fr_diff2);
metadata = vertcat(metadata1',metadata2');

clf;
figure(1);
clear g;
g=gramm('x',metadata,'y',fr_diff);

% g.stat_boxplot('notch','true')
% g.stat_violin('normalization','width');
g.stat_summary('geom',{'edge_bar','black_errorbar'},'type','ci','setylim','true');
% g.geom_jitter('width',0.6,'height',0,'dodge',0.5,'alpha',0.1);
g.set_color_options('chroma',0);
g.set_title('Change in mean FR(hz) during decoherence');
g.set_names('x','','y', 'FR Change (Hz)');
% g.axe_property('YLim',[0 8]);
g.draw();

%%
% clear g;
% figure(2);
% fr_diff = dch_avg_cellFR - before_dch_avg_cellFR;
% num_FR_increased_during_Dch = sum(fr_diff>0);
% num_FR_decreased_during_Dch = sum(fr_diff<0);
% 
% g=gramm('x',{'FR increased during decoherence','FR decreased during decoherence'},'y',[num_FR_increased_during_Dch,num_FR_decreased_during_Dch]);
% g.stat_summary('geom',{'edge_bar','black_errorbar'},'type','ci','setylim','true');
% g.set_color_options('chroma',0);
% g.axe_property('YLim',[0 2300]);
% g.set_title('Decoherence FR Increase or Decrease','FontSize',20);
% g.set_names('x',[],'y', 'Number of Cells');
% g.draw()

