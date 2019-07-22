function [all_fr, avg_all_fr, all_corrmatrix, avg_all_corrmatrix, all_corrblock, avg_all_corrblock]...
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
%% Load .mat and params files
% load specific components from .mat file
load(fullfile(matPath), 'post','posx','sp','trial'); 

% load params file

if ~exist('paramsPath','var')
    addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
    params = readtable('UniversalParams.xlsx');
else
    params = readtable(paramsPath);
end

% calculate firing rates for only good cells
cells_to_plot = sp.cids(sp.cgs==2); 


%% Preallocate and store data for all cell's firing rates
nCells = size(cells_to_plot, 2);
spatialBins = trackLength/params.SpatialBin; 
all_fr = nan(nCells, max(trial),spatialBins); % preallocate matrix of cells' firing rates across spatial bins

trials_per_block = 25;

all_corrmatrix = nan(nCells, max(trial), max(trial)); % preallocate matrix of cells' correlations between trials
all_corrblock = nan(nCells, max(trial)/trials_per_block, max(trial)/trials_per_block); % preallocate matrix of cells' correlations every 50 trials
fprintf('Calculating firing rate for %d cells with a spatial bin size of %dcm\n',nCells,params.SpatialBin);

%% Calculate firing rates
% specify inputs into the firing rates calculation
trackEnd = trackLength;
p = params;

% calculate the firing rate for a single cell across all trials

for k = 1:nCells
    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

    % get spike times and index into post for cell k 
    spike_t = sp.st(sp.clu==cells_to_plot(k));
    
    [~,~,spike_idx] = histcounts(spike_t,post);
    [~,~,spikeTrial_idx] = histcounts(spike_t,trial); % NOT SURE IF WE STILL NEED THIS
    

    % for cell k, iteratively calculate the firing rate for each trial
    for i = 1:max(trial)
        itrial_kfr = calcSmoothedFR_SpatialBin(spike_idx(i==trial(spike_idx)), posx, p, trackEnd);
        singleCellallTrialsFR(i,:) = itrial_kfr;
%         figure(1)
%         plot(itrial_kfr); %plot FR for single trial in single cell
%         figure(2) % plot rasterplot for single trial in single cell
%         plot(posx(spike_idx(i==trial(spike_idx))),trial(spike_idx(i==trial(spike_idx))),'k.');      
    end
    
%   kfr = calculateSmoothedFiringRate(spike_idx, posx, p, trackEnd); % get average firing rate collapsed across all trials
%   plot(kfr); % plot average firing rate collapsed across all trials

    % calculate trial by trial correlations for one cell 
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
    block_fr = nan(max(trial)/trials_per_block, spatialBins); % preallocate matrix
    
    % fill in preallocated matrix
    for i = 1:max(trial)/trials_per_block
        
        block_fr(i, :) = mean(singleCellallTrialsFR(trials_per_block*i-(trials_per_block-1):trials_per_block*i, :));
      
    end
    
    % calculate correlations for this block matrix
    
    corrBlock = corr(block_fr');
%     imagesc(corrBlock);
%     colorbar;
%     set(gca,'XTick',0:10:400);
%     xticklabels(xticks-100)
%     set(gca,'YTick',0:10:400);
%     yticklabels(yticks-100);
%     pause(0.5);
    
    % store one cell's firing rates and trial by trial correlation matrix
    % in a session matrix containing all cells
    all_fr(k, :, :) = singleCellallTrialsFR;
    all_corrmatrix(k, :, :) = corrMatrix;
    all_corrblock(k, :, :) = corrBlock;
    
end

avg_all_fr = squeeze(mean(all_fr, 1, 'omitnan'));
% plot(mean(avg_all_fr,2));

avg_all_corrmatrix = squeeze(mean(all_corrmatrix, 1, 'omitnan'));

avg_all_corrblock = squeeze(mean(all_corrblock, 1, 'omitnan'));

% figure(1);
% imagesc(avg_all_corrblock);
% colormap('default');
% colorbar;
% set(gca,'XTick',0:10:400);
% xticklabels(xticks-100)
% set(gca,'YTick',0:10:400);
% yticklabels(yticks-100)

% figure(1);
% rows = ceil(sqrt(nCells));
% %colormap('default');
% %colorbar;
% for i = 1:nCells
%     subplot(rows, rows, i)
%     imagesc(squeeze(all_corrblock(i, :, :)));
%     
%     set(gca,'visible','off')
%     set(gca,'XTick',[], 'YTick', [])
% end
% % 
% figure(2);
% 
% for i = 1:nCells
%     hold on
%     plot(squeeze(all_corrblock(i, 5, :)));
% end
% plot(squeeze(avg_all_corrblock(2,:)));

end