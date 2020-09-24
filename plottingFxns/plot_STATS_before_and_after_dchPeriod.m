function plot_STATS_before_and_after_dchPeriod(cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats comparing Firing Rate over during decoherence period and during equivalent length of time before dch period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
dsfactor = 100;
dch_avg_cellFR = nan(numel(cells.FRtime),1);
before_dch_avg_cellFR = nan(numel(cells.FRtime),1);

for i = 1:numel(cells.FRtime)

    fr_time = cells.FRtime(i).FRtime;
    smoothedCellFR = smoothdata(fr_time, 'sgolay',50);
    downsample_fr_time = downsample(smoothedCellFR,dsfactor);
    dch = cells.dch(i).dch;
    decoherenceTimeIdx_correctForFs = dch.decoherenceTimeIdx*dch.Fs;
    
    if ~isempty(decoherenceTimeIdx_correctForFs)
        dch_fr_time = downsample_fr_time(decoherenceTimeIdx_correctForFs); %decoherenceTimeIdx is literal time so need to multiply by Fs to get indx
        dch_avg_cellFR(i) = mean(dch_fr_time);

        dchStartIndx = decoherenceTimeIdx_correctForFs(1);
        dchIndxLength = numel(decoherenceTimeIdx_correctForFs);
        beforeStartIndx = dchStartIndx-dchIndxLength+1;
%         beforeStartIndx = 1:dchIndxLength;
        if beforeStartIndx < 1
            beforeStartIndx = 1;
        end
        before_dch_fr_time = downsample_fr_time(beforeStartIndx:dchStartIndx);
        before_dch_avg_cellFR(i) = mean(before_dch_fr_time);
    else
        fprintf('No Decoherence period\n') 
    end

end

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
g=gramm('x',metadata,'y',fr,'color',color);
% g.stat_violin('normalization','width','dodge',0,'fill','edge','half','true');
g.stat_boxplot('notch','true');
g.set_title('FR(hz) Control vs Decoherence Period');
g.set_names('x','Mice','y', 'Hz');
g.draw()

%%
close all; clear g;
fr_diff = abs(dch_avg_cellFR - before_dch_avg_cellFR);
nanIndx = isnan(fr_diff);
fr_diff = fr_diff(~nanIndx);

metadata = cells.metadata(:,2);
metadata = metadata(~nanIndx);

g=gramm('x',metadata,'y',fr_diff);
% g.geom_jitter('width',0.6,'height',0,'dodge',1);
% g.stat_boxplot('notch','true')
g.stat_summary('geom',{'edge_bar','black_errorbar'},'type','ci','setylim','true')
g.set_title('Change in mean FR(hz) during decoherence');
g.set_names('x','Mice','y', 'Abs(Change Hz)');
g.draw();

%%


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