function [post,posx,sp, all_fr, avg_all_fr, all_corrmatrix, avg_all_corrmatrix,...
    all_waveforms, cells_to_plot,spike_depth,...
    all_correlationScore, trial,all_cellCorrScore,trials_corrTemplate, avg_all_cellCorrScore, avg_cell_fr]...
    = calcFRmapCorrMatrixAllCells(matPath, trackLength, paramsPath)

% John Wen 7/1/19
% Kei Masuda 7/3/19
% Runs Malcolm's firing rate code on all good cells within a .mat file. 
% Calculate firing rate by dividing spike counts by occupancy
% inputs:
%     matPath: path to .mat file after running sync_vr_to_np.m. Specify as
%              string.
%              e.g. '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190620_keicontrasttrack_ketamine1_g0/G4_190620_keicontrasttrack_baseline+cntrlinjx+ketamine'
%     trackLength: specify the length of the track
%     paramsPath: Optional argument. Path to parameters file, which includes 
%                 spatial bin size and smoothing kernel. 
%                 Default = '/Volumes/groups/giocomo/export/data/Projects/ ...
%                 JohnKei_NPH3/UniversalParams'
% outputs:
%     all_fr: smoothed firing rate over position for each cell 
%             (cells x trials x spatial bins)
%     avg_all_fr: smoothed firing rate over position, averaged over cells
%                 (trials x spatial bins)
%     all_corrmatrix: trial by trial correlation for all cells (cells x
%     trial x trial)
%     avg_all_corrmatrix: mean of all_corrmatrix across cells
%     all_corrblock: correlation across every block of 50 trials (cells x
%     trial/50 x trial/50)
%     avg_all_corrblock: average of all_corrblock across cells
%%
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));

% addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
% addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));
%% Load .mat and params files
% load specific components from .mat file
load(fullfile(matPath), 'post','posx','sp','trial'); 

if ~exist('paramsPath','var')
    addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
    params = readtable('UniversalParams.xlsx');
else
    params = readtable(paramsPath);
end

trackEnd = trackLength;
p = params;

% calculate firing rates for only good cells
cells_to_plot = sp.cids(sp.cgs==2); 
% compute some useful information (like spike depths)
[~, spikeDepths, ~, ~, ~, ~, waveforms] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);


% get spike depths
spike_depth = nan(numel(cells_to_plot),1);
for k = 1:numel(cells_to_plot)
    spike_depth(k) = median(spikeDepths(sp.clu==cells_to_plot(k)));
end


% sort cells_to_plot by spike_depth (descending)
[spike_depth,sort_idx] = sort(spike_depth,'descend');
cells_to_plot = cells_to_plot(sort_idx);

%% Calculate speed and filter out stationary periods (speed<2cm)
% speed = calcSpeed(posx, p);
% [trial,posx,post] = speedFilterData(trial,posx, post, speed);

%% Preallocate and store data for all cell's firing rates
nCells = size(cells_to_plot, 2);
spatialBins = trackLength/params.SpatialBin; 
all_fr = nan(nCells, max(trial),spatialBins); % preallocate matrix of cells' firing rates across spatial bins

trials_per_block = 25;
trials_corrTemplate = 50;

all_corrmatrix = nan(nCells, max(trial), max(trial)); % preallocate matrix of cells' correlations between trials
% all_corrblock = nan(nCells, floor(max(trial)/trials_per_block), floor(max(trial)/trials_per_block)); % preallocate matrix of cells' correlations every 50 trials
all_waveforms= nan(nCells, size(waveforms,2)); 
all_cellCorrScore = nan(nCells, numel(trials_corrTemplate:max(trial)));
all_correlationScore = nan(nCells, 2);
fprintf('Calculating firing rate for %d cells with a spatial bin size of %dcm\n',nCells,params.SpatialBin);

%% Calculate firing rates
% specify inputs into the firing rates calculation

for k = 1:nCells
    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));
    
    % get spike times and index into post for cell k 
    spike_t = sp.st(sp.clu==cells_to_plot(k));
    
    [~,~,spike_idx] = histcounts(spike_t,post);
    spike_idx(spike_idx==0) = []; %remove spike indexes that don't exist in post (e.g. spikes that happen while the animal is stationary when post was speed filtered)
    
    % for cell k, iteratively calculate the firing rate for each trial
%     figure(5);% plot rasterplot for single trial in single cell
%     clf;
    singleCellallTrialsFR = nan(max(trial),spatialBins);
    for i = 1:max(trial)
        itrial_kfr = calcSmoothedFR_SpatialBin(spike_idx(i==trial(spike_idx)), posx(i==trial),posx, p, trackEnd);
%         itrial_kfr = calcSmoothedFR_SpatialBin_speedFiltered(spike_t, post,posx, p, trackEnd,speed, trial, i);
        singleCellallTrialsFR(i,:) = itrial_kfr;
%         figure(1)
%         plot(itrial_kfr); %plot FR for single trial in single cell
%         plot(posx(spike_idx(i==trial(spike_idx))),trial(spike_idx(i==trial(spike_idx))),'k.');hold on;
    end
    
    
%   kfr = calculateSmoothedFiringRate(spike_idx, posx, p, trackEnd); % get average firing rate collapsed across all trials
%   plot(kfr); % plot average firing rate collapsed across all trials

    % for cell k, calculate correlation template from 1st 40trials
    corrTemplate = mean(singleCellallTrialsFR(1:trials_corrTemplate,:));
%     plot(corrTemplate)
    cellCorrScore = nan(numel(trials_corrTemplate:max(trial)),1);
    
    for i = trials_corrTemplate:max(trial)
        corrMatrix = corrcoef(corrTemplate, singleCellallTrialsFR(i,:));
        corrScore = corrMatrix(1,2);
        cellCorrScore(i-trials_corrTemplate+1) = corrScore;
    end
    
    cellCorrScore = fillmissing(cellCorrScore,'spline');
    smoothCellCorrScore= sgolayfilt(cellCorrScore, 3, 51);
    corrPeak_preKet = max(findpeaks(smoothCellCorrScore(25:50)));
    corrTrough_postKet = -max(findpeaks(-smoothCellCorrScore(50:75)));
    correlationScore = corrPeak_preKet - corrTrough_postKet;
    
%     plot(cellCorrScore)
    all_cellCorrScore(k,:) = cellCorrScore;
    all_correlationScore(k,:) = [correlationScore, spike_depth(k)];

    %% calculate trial by trial correlation matrix for one cell 
    corrMatrix = corr(singleCellallTrialsFR');
    
    % CODE to plot correlation matrix per trial
    %figure(1);
%     imagesc(corrMatrix);
%     colorbar;
%     set(gca,'XTick',0:10:400);
%     xticklabels(xticks-100)
%     set(gca,'YTick',0:10:400);
%     yticklabels(yticks-100)

    % calculate correlations across blocks of trials (every 50 trials)
    
    % create trial-averaged firing rate matrix (average every 50 trials)
    numBlocks = ceil(max(trial)/trials_per_block);
    block_fr = nan(numBlocks, spatialBins); % preallocate matrix
    
%     % fill in preallocated matrix
%     for i = 1:numBlocks
%         try
%             block_fr(i, :) = mean(singleCellallTrialsFR(trials_per_block*i-(trials_per_block-1):trials_per_block*i, :));
%         catch
%             block_fr(i, :) = NaN;
%         end
%     end
    
    % calculate correlations for this block matrix
%     block_fr(any(isnan(block_fr), 2), :) = [];
%     corrBlock = corr(block_fr');
%     imagesc(corrBlock);
%     colorbar;
%     set(gca,'XTick',0:numBlocks:max(trial));
%     xticklabels(xticks-100)
%     set(gca,'YTick',0:numBlocks:max(trial));
%     yticklabels(yticks-100);
%     pause(0.5);
%     
    % store one cell's firing rates and trial by trial correlation matrix
    % in a session matrix containing all cells
    all_fr(k, :, :) = singleCellallTrialsFR;
    all_corrmatrix(k, :, :) = corrMatrix;
%     all_corrblock(k, :, :) = corrBlock;
    all_waveforms(k,:) = waveforms(sp.cids==cells_to_plot(k),:);

end
%%
avg_all_fr = squeeze(mean(all_fr, 1, 'omitnan'));
% plot(mean(avg_all_fr,2));
avg_cell_fr = squeeze(mean(mean(all_fr, 2, 'omitnan'),3));

avg_all_corrmatrix = squeeze(mean(all_corrmatrix, 1, 'omitnan'));

% avg_all_corrblock = squeeze(mean(all_corrblock, 1, 'omitnan'));


avg_all_cellCorrScore = squeeze(mean(all_cellCorrScore, 1, 'omitnan'));

% PLOT AVERAGE CORRELATION BLOCK MATRIX
% figure(1);
% imagesc(avg_all_corrblock);
% colormap('default');
% colorbar;
% set(gca,'XTick',0:25:400);
% xticklabels(xticks-100)
% set(gca,'YTick',0:25:400);
% yticklabels(yticks-100)
% 
% PLOT CORRELATION BLOCK MATRIX FOR EACH CELL
% figure(2);
% set(gcf,'Position',[100 100 1000 1000]); 
% rows = ceil(sqrt(nCells));
% %colormap('default');
% %colorbar;
% for i = 1:nCells
%     subplot(rows, rows, i)
%     imagesc(squeeze(all_corrblock(i, :, :)));
% %     title(sprintf('c%d, d=%d',cells_to_plot(i),round(spike_depth(i))));
%     set(gca,'visible','off')
%     set(gca,'XTick',[], 'YTick', [])
% end
% 
% PLOT WAVEFORMS FOR EACH CELL
% figure(3);
% set(gcf,'Position',[100 100 4000 1000]); 
% rows = ceil(sqrt(nCells));
% for i = 1:nCells
%     subplot(rows, rows, i)
%     plot(all_waveforms(i, :));
%     title(sprintf('c%d, d=%d',cells_to_plot(i),round(spike_depth(i))));
% %     set(gca,'visible','off')
%     set(gca,'XTick',[], 'YTick', [])
% end


% PLOT SINGLE ROW OF EACH CORR MATRIX
% figure(3);
% for i = 1:nCells
%     hold on
%     plot(squeeze(all_corrblock(i, 5, :)));
% end
% plot(squeeze(avg_all_corrblock(2,:)));



% PLOT FR for each cell
% figure(5);
% set(gcf,'Position',[100 100 1000 1000]); 
% rows = ceil(sqrt(nCells));
% %colormap('default');
% %colorbar;
% for i = 1:nCells
%     subplot(rows, rows, i)
%     clear g;
%     x = -99:200;
%     y = squeeze(all_fr(i, :,:));
%     plot(x,mean(y'));
%     ttl = title(sprintf('c%d, d=%d',cells_to_plot(i),round(spike_depth(i))));
%     ttl.FontSize = 10;
%     set(gca,'XTick',[], 'YTick', [])
% end
% 

%%
% %% PLOT CORRLEATION PLOT FOR ONE CELL
% figure(6);
% hold on;
% p1 = plot((trials_corrTemplate:max(trial)),cellCorrScore, 'b', 'LineWidth',0.05);
% p1.Color(4) = 0.15;
% smoothY= sgolayfilt(cellCorrScore, 3, 25);
% plot((trials_corrTemplate:max(trial)), smoothY,'r','LineWidth',3); 
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('G4-190620: Cell 85 Correlation Score')
% xlabel('Trials')
% ylabel('Correlation Score compared to First 50 Trials')
% 
% %% PLOT CORRLEATION PLOT FOR SESSION
% figure(7);
% hold on;
% p1 = plot((trials_corrTemplate:max(trial)),avg_all_cellCorrScore, 'b', 'LineWidth',0.05);
% p1.Color(4) = 0.15;
% smoothY= sgolayfilt(avg_all_cellCorrScore, 3, 25);
% plot((trials_corrTemplate:max(trial)), smoothY,'r','LineWidth',3); 
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('G4-190620: Average Correlation Score')
% xlabel('Trials')
% ylabel('Correlation Score compared to First 50 Trials')
% 
% %% PLOT Distribution of CORR SCORE for all cells
% figure(8);
% hold on;
% p1 = plot((trials_corrTemplate:max(trial)),avg_all_cellCorrScore, 'b', 'LineWidth',0.05);
% p1.Color(4) = 0.15;
% smoothY= sgolayfilt(avg_all_cellCorrScore, 3, 25);
% plot((trials_corrTemplate:max(trial)), smoothY,'r','LineWidth',3); 
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('G4190620 - Avg Correlation Score')
% xlabel('Trials')
% ylabel('Correlation Score compared to First 49 Trials')
% 
% %% PLOT CORRELATION TEMPLATE FOR ONE CELL
% figure(9);
% hold on;
% plot(corrTemplate,'r','LineWidth',2); 
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('Baseline Correlation Template')
% xlabel('Position Bin (Each is 2cm)')
% ylabel('Avg Baseline Firing Rate (HZ)')
% 
% %%
% %% PLOT CORRELATION LINE FOR ONE CELL
% figure(10);
% hold on;
% plot(cellCorrScore,'r','LineWidth',2); 
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('Correlation to Baseline Template')
% xlabel('Trials')
% ylabel('Correlation to Baseline Template')
% 
% %% Plot Distribution of Correlation Scores
% figure(11);
% hold on;
% histfit(all_correlationScore(:,1),50, 'kernel')
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('Distribution of Ketamine Correlation Effect Scores')
% xlabel('Pre vs Post Ketamine Correlation Effect Score')
% ylabel('Number of Cells')
% 
% %% PLOT SCATTER of Corr Scorr vs Cell Depth
% figure(12);
% hold on;
% scatter(all_correlationScore(:,2),all_correlationScore(:,1), 150,'filled','r');
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('Ketamine Correlation Effect Score vs Cell Depth')
% xlabel('Distance from Tip of Probe')
% ylabel('Pre vs Post Ketamine Correlation Effect Score')
% 
% %% PLOT Scatter of Ketamien Corr Score vs FR
% figure(12);
% hold on;
% scatter(avg_cell_fr,all_correlationScore(:,1), 150,'filled','r');
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('Ketamine Correlation Effect Score vs Firing Rate')
% xlabel('FR of Cell (Hz)')
% ylabel('Pre vs Post Ketamine Correlation Effect Score')
end