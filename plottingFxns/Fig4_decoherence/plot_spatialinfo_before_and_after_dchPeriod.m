function plot_spatialinfo_before_and_after_dchPeriod(cells)

% function constants 
plotWidth = 500;
plotHeight = 1000;
spatialBinNum = 200;
MPP = 0.125; %Minimum Peak Prominence

numCells = size(cells.metadata,1);
drawFig = true;


% declare global vaariables
all_Ispike_cntrl = NaN(numCells,1);
all_Ispike_dch = NaN(numCells,1);

if drawFig
    figure('Position',[200 200 plotWidth plotHeight]); hold on;
end
%%
for i = 1:numCells

    
    singleCellSpatialFR = squeeze(cells.spatialFR2(i,:,:));
    

    % Get individual cells' posx,post,spatialBin, speed values for spatial
    % info calcuation
    dch = cells.dch(i).dch;
    posx = cells.posX(i).posx;
    post = cells.posT(i).post;
    speed = cells.speed(i).speed;
    trial = cells.trial(i).trial;
    spatialBin = 2; %2cm
    decoherenceTrialIdx = dch.decoherenceIdx;
    
    if ~isempty(decoherenceTrialIdx)
        try 
            dchIndxLength = numel(decoherenceTrialIdx);
            controlTrials = 101-dchIndxLength:100;
%             controlTrials = 1:dchIndxLength;
%             if controlTrials(1) < 1
%                 controlTrials = 1:dchIndxLength;
%             end
%             control_spatialFR = singleCellSpatialFR(controlTrials);
%             dch_spatialFR = singleCellSpatialFR(decoherenceTrialIdx);
                   
    
            %% Calculate Spatial Information for dch and cntrl
        
            trialBlockSpatialInformation = nan(max(trial),2);
            for j = 1:max(trial)
                trial_posx = posx(trial==j);
                trial_post = post(trial==j);
                trial_speed = speed(trial==j);
                linearFractionalOccupancyBlock = calculate_1D_LFO(trial_posx,trial_post,spatialBin,trial_speed);%calculates linear fraction occupancy for a trial 
        
                [I_sec_block, I_spike_block] = calculate_1DspatialInformation(singleCellSpatialFR(j,:),linearFractionalOccupancyBlock);
                trialBlockSpatialInformation(j,:) = [I_sec_block, I_spike_block];
            end

            
            I_spike_cntrl = nanmean(trialBlockSpatialInformation(controlTrials,1));
            I_spike_dch = nanmean(trialBlockSpatialInformation(decoherenceTrialIdx,1));
            %% Save values to global variable
            all_Ispike_cntrl(i) = I_spike_cntrl;
            all_Ispike_dch(i) = I_spike_dch;
            fprintf('finished %i iterations\n', i);
        catch
            fprintf('FAILED #%i iteration\n', i)
        end
    else
        fprintf('No Decoherence period\n') 
    end
end


%% Summary Plots

[h,p] = ttest2(all_Ispike_cntrl, all_Ispike_dch);
fprintf('The hypothesis test is %i; P-value is %.02d\n',h,p);
 
% Compare Before Dch and After Dch Figures
close all; 
data = vertcat(all_Ispike_cntrl,all_Ispike_dch);
nanIndx = isnan(data);
data = data(~nanIndx);

% % Get Metadata Data
% metadata = cells.metadata(:,2);
% metadata = vertcat(metadata,metadata);
% metadata = metadata(~nanIndx);

color = vertcat(...
    repmat({'Control'}, numCells),...
    repmat({'Decoherence'}, numCells)...
    );
color = color(~nanIndx);
%
clear g;
% g=gramm('x',data,'y',color);
g=gramm('x',color,'y',data,'color',color);
% g.stat_bin();
g.stat_violin('normalization','width')

% g.stat_boxplot();
% g.stat_boxplot('notch','true','geom_jitter');
% g.stat_summary();
% g.stat_bin('geom','stacked_bar')
g.set_title('Spatial Information Score Control vs Decoherence Period');
g.set_names('x','','y', 'SI Score');
customColorMap = [ % 0.5 0.5 0.5 %grey
    0.8 0.2 0.8 %magenta
    0 0.8 0.2]; %green
g.set_color_options('map',customColorMap);
% g.geom_jitter();
g.draw()

%%
% Difference in SI before Dch and after Dch across mice
% close all; 
% data = all_Ispike_cntrl-all_Ispike_dch;
% nanIndx = isnan(data);
% data = data(~nanIndx);
% 
% % Get Mice Name Data
% metadata = cells.metadata(:,2);
% metadata = vertcat(metadata,metadata);
% metadata = metadata(~nanIndx);
% 
% clear g;
% g=gramm('x',metadata,'y',data);
% g.stat_bin('geom','stacked_bar')
% % g.stat_bin();
% % g.stat_summary();
% % g.geom_jitter();
% % g.stat_boxplot('notch','true','geom_jitter');
% % g.stat_summary();
% g.set_title('Spatial Information Score Control vs Decoherence Period');
% g.set_names('x','','y', 'Count');
% customColorMap = [ 0.5 0.5 0.5 %grey
%     0.8 0.2 0.8 %magenta
%     0 0.8 0.2]; %green
% g.set_color_options('map',customColorMap);
% g.draw()