function plot_cells1_vs_cells2_dchPeriod_STATS(cells1,cells2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats comparing Firing Rate over during decoherence period and during equivalent length of time before dch period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
dchAndBefore1 = find_dch_and_equivalent_beforetime(cells1);
dchAndBefore2 = find_dch_and_equivalent_beforetime(cells2);
%%
abs_fr_diff1 = abs(dchAndBefore1.dch_avg_cellFR - dchAndBefore1.before_dch_avg_cellFR);
abs_fr_diff2 = abs(dchAndBefore2.dch_avg_cellFR - dchAndBefore2.before_dch_avg_cellFR);

clear g;
figure();


nanIndx1 = isnan(abs_fr_diff1);
abs_fr_diff1 = abs_fr_diff1(~nanIndx1);
metadata1 = cells1.metadata(:,4);
metadata1 = metadata1(~nanIndx1);


nanIndx2 = isnan(abs_fr_diff2);
abs_fr_diff2 = abs_fr_diff2(~nanIndx2);
metadata2 = cells2.metadata(:,4);
metadata2 = metadata2(~nanIndx2);


[h,p] = ttest2(abs_fr_diff1, abs_fr_diff2);
fprintf('The hypothesis test is %i; P-value is %.02d\n',h,p);


abs_fr_diff = horzcat(abs_fr_diff1,abs_fr_diff2);
metadata = vertcat(metadata1,metadata2);


g=gramm('x',metadata,'y',abs_fr_diff);
% g.geom_jitter('width',0.6,'height',0,'dodge',1);
% g.stat_boxplot('notch','true')
% g.stat_violin('normalization','width')
g.stat_summary('geom',{'edge_bar','black_errorbar'},'type','ci','setylim','true');
g.set_color_options('chroma',0);
g.set_title('Absolute change in mean FR(hz) during decoherence');
g.set_names('x','Mice','y', 'Abs(Change Hz)');
g.axe_property('YLim',[0 8]);
g.draw();

%%
clear g;
figure();
fr_diff = dch_avg_cellFR - before_dch_avg_cellFR;
num_FR_increased_during_Dch = sum(fr_diff>0);
num_FR_decreased_during_Dch = sum(fr_diff<0);

g=gramm('x',{'FR increased during decoherence','FR decreased during decoherence'},'y',[num_FR_increased_during_Dch,num_FR_decreased_during_Dch]);
g.stat_summary('geom',{'edge_bar','black_errorbar'},'type','ci','setylim','true');
g.set_color_options('chroma',0);
g.axe_property('YLim',[0 2300]);
g.set_title('Decoherence FR Increase or Decrease','FontSize',20);
g.set_names('x',[],'y', 'Number of Cells');
g.draw()

