function plot_spatialinfo_before_and_after_dchPeriod(cells)

% function constants 
plotWidth = 500;
plotHeight = 1000;
spatialBinNum = 200;
MPP = 0.125; %Minimum Peak Prominence

numCells = size(cells.metadata,1);
drawFig = true;


% declare global vaariables
all_Isec_cntrl = NaN(numCells,1);
all_Isec_dch = NaN(numCells,1);

if drawFig
    figure('Position',[200 200 plotWidth plotHeight]); hold on;
end
%%
for i = 1:numCells

    smoothFR = cells.spatialFRsmooth;
    spatialInfo = cells.spatialInfo;

    % Get individual cells' posx,post,spatialBin, speed values for spatial
    % info calcuation
    dch = cells.dch(i).dch;
    posx = cells.posX(i).posx;
    post = cells.posT(i).post;
    spatialBin = 2; %2cm
    decoherenceTrialIdx = dch.decoherenceIdx;
    
    if ~isempty(decoherenceTrialIdx)
        try 
            dchIndxLength = numel(decoherenceTrialIdx);
    %         controlTrials = 101-dchIndxLength:100;
            controlTrials = 1:dchIndxLength;
            if controlTrials(1) < 1
                controlTrials = 1:dchIndxLength;
            end
            control_spatialFR = squeeze(smoothFR(i,controlTrials,:));
            dch_spatialFR = squeeze(smoothFR(i,decoherenceTrialIdx,:));
                   
    
            %% Calculate Spatial Information for dch and cntrl
            linearFractionalOccupancy = calculate_1D_LFO(posx,post,spatialBin, speed);
            [I_sec_cntrl, I_spike_cntrl] = calculate_1DspatialInformation(control_spatialFR,linearFractionalOccupancy);
            [I_sec_dch, I_spike_dch] = calculate_1DspatialInformation(dch_spatialFR,linearFractionalOccupancy);
    
    
            %% Save values to global variable
            all_Isec_cntrl(i) = I_sec_cntrl;
            all_Isec_dch(i) = I_sec_dch;
            fprintf('finished %i iterations\n', i);
        catch
            fprintf('FAILED #%i iteration\n', i)
        end
    else
        fprintf('No Decoherence period\n') 
    end
end


%% Summary Plots

[h,p] = ttest2(all_Isec_cntrl, all_Isec_dch);
fprintf('The hypothesis test is %i; P-value is %.02d\n',h,p);
% Compare Before Dch and After Dch Figures
close all; 
data = vertcat(all_Isec_cntrl,all_Isec_dch);
nanIndx = isnan(data);
data = data(~nanIndx);

% Get Mice Name Data
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
% g=gramm('x',metadata,'y',data,'color',color,'subset',metadata);
g=gramm('x',metadata);
% g.stat_bin();
% g.stat_violin();
% g.stat_boxplot();
% g.stat_boxplot('notch','true','geom_jitter');
% g.stat_summary();
g.stat_bin('geom','stacked_bar')
g.set_title('Spatial Information Score Control vs Decoherence Period');
g.set_names('x','','y', 'Hz');
customColorMap = [ % 0.5 0.5 0.5 %grey
    0.8 0.2 0.8 %magenta
    0 0.8 0.2]; %green
g.set_color_options('map',customColorMap);
g.draw()

%%
% Difference in SI before Dch and after Dch across mice
close all; 
data = all_Isec_cntrl-all_Isec_dch;
nanIndx = isnan(data);
data = data(~nanIndx);

% Get Mice Name Data
metadata = cells.metadata(:,2);
metadata = vertcat(metadata,metadata);
metadata = metadata(~nanIndx);

clear g;
g=gramm('x',metadata,'y',data);
g.stat_bin('geom','stacked_bar')
% g.stat_bin();
% g.stat_summary();
% g.geom_jitter();
% g.stat_boxplot('notch','true','geom_jitter');
% g.stat_summary();
g.set_title('Spatial Information Score Control vs Decoherence Period');
g.set_names('x','','y', 'Count');
customColorMap = [ 0.5 0.5 0.5 %grey
    0.8 0.2 0.8 %magenta
    0 0.8 0.2]; %green
g.set_color_options('map',customColorMap);
g.draw()