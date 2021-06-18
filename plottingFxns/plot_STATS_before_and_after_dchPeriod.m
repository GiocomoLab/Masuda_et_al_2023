function plot_STATS_before_and_after_dchPeriod(cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats comparing Firing Rate over during decoherence period and during equivalent length of time before dch period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
% dsfactor = 100;
% 
% % initialize variables
% dch_avg_cellFR = nan(size(cells.FRtime));
% before_dch_avg_cellFR = nan(size(cells.FRtime));
% 
% for i = 1:numel(cells.FRtime)
% 
%     fr_time = cells.FRtime(i).FRtime;
%     smoothedCellFR = smoothdata(fr_time, 'gaussian',50);
%     downsample_fr_time = downsample(smoothedCellFR,dsfactor);
%     dch = cells.dch(i).dch;
%     decoherenceTimeIdx_correctForFs = dch.decoherenceTimeIdx*dch.Fs; %decoherenceTimeIdx is literal time so need to multiply by Fs to get indx
%     
%     if ~isempty(decoherenceTimeIdx_correctForFs)
%         dch_fr_time = downsample_fr_time(decoherenceTimeIdx_correctForFs); 
%         dch_avg_cellFR(i) = mean(dch_fr_time);
% 
%         dchStartIndx = decoherenceTimeIdx_correctForFs(1);
%         dchIndxLength = numel(decoherenceTimeIdx_correctForFs);
%         beforeStartIndx = dchStartIndx-dchIndxLength+1;
% %         beforeStartIndx = 1:dchIndxLength;
%         if beforeStartIndx < 1
%             beforeStartIndx = 1;
%         end
%         before_dch_fr_time = downsample_fr_time(beforeStartIndx:dchStartIndx);
%         before_dch_avg_cellFR(i) = mean(before_dch_fr_time);
%     end
% 
% end
dchAndBefore = find_dch_and_equivalent_beforetime(cells);
before_dch_avg_cellFR = dchAndBefore.before_dch_avg_cellFR;
dch_avg_cellFR = dchAndBefore.dch_avg_cellFR;

%%
[h,p] = ttest2(before_dch_avg_cellFR, dch_avg_cellFR);
fprintf('The hypothesis test is %i; P-value is %.02d\n',h,p);
%% Compare Before and After Figures
close all; 
fr = vertcat(before_dch_avg_cellFR,dch_avg_cellFR);
nanIndx = isnan(fr);
fr = fr(~nanIndx);

metadata = cells.metadata(:,2);
metadata = vertcat(metadata,metadata);
metadata = metadata(~nanIndx);

color = vertcat(...
    repmat({'Control'}, size(before_dch_avg_cellFR)),...
    repmat({'Decoherence'}, size(dch_avg_cellFR))...
    );
color = color(~nanIndx);
%
clear g;
g=gramm('x',color,'y',fr,'color',color);
% g.stat_violin('normalization','width','dodge',0,'fill','edge','half','true');
g.stat_boxplot('notch','true');
g.set_title('FR(hz) Control vs Decoherence Period');
g.set_names('x','Mice','y', 'Hz');
g.draw()

%%
clear g;
figure();
abs_fr_diff = abs(dch_avg_cellFR - before_dch_avg_cellFR);
nanIndx = isnan(abs_fr_diff);
abs_fr_diff = abs_fr_diff(~nanIndx);

metadata = cells.metadata(:,2);
metadata = metadata(~nanIndx);

g=gramm('x',metadata,'y',abs_fr_diff);
% g.geom_jitter('width',0.6,'height',0,'dodge',1);
% g.stat_boxplot('notch','true')
% g.stat_violin('normalization','width')
g.stat_summary('geom',{'edge_bar','black_errorbar'},'type','ci','setylim','true')
g.set_color_options('chroma',0)
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

