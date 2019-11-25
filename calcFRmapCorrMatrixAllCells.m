function [post,posx,sp, all_fr, avg_all_fr, all_corrmatrix, avg_all_corrmatrix,...
    all_waveforms, cells_to_plot,spike_depth,...
    all_drugEffectScores, trial,all_cellCorrScore,trials_corrTemplate, avg_all_cellCorrScore, avg_cell_fr,...
    trial_ds, all_frTime]...
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
%     post, pox, sp: pass through vars from spikes
%     all_fr: smoothed firing rate over position tensor 
%             (cells x trials x spatial bins)
%     avg_all_fr: smoothed firing rate over position, averaged over cells
%                 (trials x spatial bins)
%     all_corrmatrix: trial by trial correlation for all cells (cells x
%     trial x trial)
%     avg_all_corrmatrix: mean of all_corrmatrix across cells
%     all_waveforms: template of all cell's waveforms
%     cells_to_plot: indices of good cells to plot
%     spiked_depth: spike depth of each cell
%     all_drugEffectScores = [drugFRdiff,cntrlFRdiff, drugFREffectScore, drugCorrEffectScore, spike_depth(k)]; 
%     trial: trial information about every sp sample
%     all_cellCorrScore: cellCorrScore curves compared to baseline template
%     trials_corrTemplate: baseline template for each cell
%     avg_all_cellCorrScore: avarge cellCorrScores curves for every session
%     avg_cell_fr: average fr for all cells across spatial bins
%     trials_ds: downsampled trial vector matching frTime
%     all_frTime: smooted firing rate over time 
%%
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));

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

% pull out good cells and sort by spike depths
[cells_to_plot, spike_depth,waveforms] = spPreProcess(sp);

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
all_cellCorrScore = nan(nCells, numel(1:max(trial)));
all_cellStabilityScore = nan(nCells, numel(1:max(trial)));
all_drugEffectScores = nan(nCells, 5); %FR score, drug correlation effect score, spikeDepth
fprintf('Calculating firing rate for %d cells with a spatial bin size of %dcm\n',nCells,params.SpatialBin);
%%

% calculate firing rate by time
all_frTime = [];
ds_factor = 1;
smoothSigma = 10;
trial_ds = downsample(trial, ds_factor); 


%% Calculate firing rates
% specify inputs into the firing rates calculation

for k = 1:nCells
    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));
    
    % get spike times and index into post for cell k 
    spike_t = sp.st(sp.clu==cells_to_plot(k));
    
    % calculate firing rate by time
    fr = calcSmoothedFR_Time(post, spike_t, ds_factor, smoothSigma);
    all_frTime(k, :) = fr;
    
     % Bin spikes into time associated spatial bins
    [~,~,spike_idx] = histcounts(spike_t,post);
    spike_idx(spike_idx==0) = []; %remove spike indexes that don't exist in post (e.g. spikes that happen while the animal is stationary when post was speed filtered)
    
    % for cell k, iteratively calculate the firing rate for each trial
    singleCellallTrialsFR = nan(max(trial),spatialBins);
    for i = 1:max(trial)
        itrial_kfr = calcSmoothedFR_SpatialBin(spike_idx(i==trial(spike_idx)), posx(i==trial),posx, p, trackEnd);
        singleCellallTrialsFR(i,:) = itrial_kfr;
    end
    %% 
    [drugCorrEffectScore, cellCorrScore, corrTemplate] = calculateCorrScore(singleCellallTrialsFR, trials_corrTemplate);    
    all_cellCorrScore(k,:) = cellCorrScore;
    
    all_cellStabilityScore(k,:) = cellStabilityScore;
    
    [drugFRdiff,cntrlFRdiff,drugFREffectScore] = calculateFRScore(singleCellallTrialsFR, 100, 50);
    
    all_drugEffectScores(k,:) = [drugFRdiff,cntrlFRdiff, drugFREffectScore, drugCorrEffectScore, spike_depth(k)];

    %% calculate trial by trial correlation matrix for one cell 
    corrMatrix = corr(singleCellallTrialsFR');
    
    % calculate correlations across blocks of trials (every 50 trials)
    
    % create trial-averaged firing rate matrix (average every 50 trials)
    numBlocks = ceil(max(trial)/trials_per_block);
    % block_fr = nan(numBlocks, spatialBins); % preallocate matrix
    
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

end