function plot_correlation_before_and_after_dchPeriod(cells)

% function constants 
plotWidth = 500;
plotHeight = 1000;
spatialBinNum = 200;
MPP = 0.125; %Minimum Peak Prominence

numCells = size(cells.metadata,1);
drawFig = true;


% declare global vaariables
all_pcntrl = NaN(numCells,1);
all_pdch = NaN(numCells,1);
all_cells_corrcoef = NaN(numCells,1);

if drawFig
    figure('Position',[200 200 plotWidth plotHeight]); hold on;
end
%%
for i = 1:numCells

    smoothFR = cells.spatialFRsmooth;
    spatialInfo = cells.spatialInfo;


    dch = cells.dch(i).dch;

    decoherenceTrialIdx = dch.decoherenceIdx;
    
    if ~isempty(decoherenceTrialIdx)
        dchIndxLength = numel(decoherenceTrialIdx);
%         controlTrials = 101-dchIndxLength:100;
        controlTrials = 1:dchIndxLength;
        if controlTrials(1) < 1
            controlTrials = 1:dchIndxLength;
        end
        control_spatialFR = squeeze(smoothFR(i,controlTrials,:));
        dch_spatialFR = squeeze(smoothFR(i,decoherenceTrialIdx,:));
               
        flat_control_spatialFR = reshape(control_spatialFR.',1,[]);
        flat_dch_spatialFR = reshape(dch_spatialFR.',1,[]); 

        flat_control_spatialFR = fillmissing(flat_control_spatialFR,'previous');
        flat_dch_spatialFR = fillmissing(flat_dch_spatialFR,'previous');
        
        %% Correlate
         [rhocoef,pval] = corrcoef(flat_control_spatialFR, flat_dch_spatialFR);
           
       
        %% Save values to global variable
        all_pcntrl(i) = p_cntrl;
        all_pdch(i) = p_dch;
        all_cells_corrcoef(i) = rhocoef(2,1);
        fprintf('finished %i iterations\n', i);   
    else
        fprintf('No Decoherence period\n') 
    end
end


%% Summary Plots

[h,p] = ttest2(all_pcntrl, all_pdch);
fprintf('The hypothesis test is %i; P-value is %.02d\n',h,p);
% Compare Before Dch and After Dch Figures
close all; 
data = vertcat(all_pcntrl,all_pdch);
nanIndx = isnan(data);
data = data(~nanIndx);

metadata = cells.metadata(:,2);
metadata = vertcat(metadata,metadata);
metadata = metadata(~nanIndx);

color = vertcat(...
    repmat({'Control'}, numCells),...
    repmat({'Decoherence'}, numCells)...
    );
color = color(~nanIndx);
%
clear g;
g=gramm('x',color,'y',data,'color',color);
% g.stat_bin();
% g.stat_violin('normalization','width','dodge',0,'fill','edge','half','true');
g.stat_boxplot('notch','true','geom_jitter');
% g.stat_summary();
g.set_title('Autocorrelation Score Control vs Decoherence Period');
g.set_names('x','','y', 'Hz');
customColorMap = [ % 0.5 0.5 0.5 %grey
    0.8 0.2 0.8 %magenta
    0 0.8 0.2]; %green
g.set_color_options('map',customColorMap);
g.draw()