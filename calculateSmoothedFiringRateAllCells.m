function all_fr = calculateSmoothedFiringRateAllCells(matPath, trackLength, paramsPath)

% John Wen 7/1/19
% Kei Masuda 7/3/19
% Runs Malcolm's firing rate code on all good cells within a .mat file. 
% Calculate firing rate by dividing spike counts by occupancy

% inputs:
%     matPath: path to .mat file after running sync_vr_to_np.m. Specify as
%             string.
%       e.g. '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190620_keicontrasttrack_ketamine1_g0/G4_190620_keicontrasttrack_baseline+cntrlinjx+ketamine'
%     trackLength: specify the length of the track
%     paramsPath: Optional argument. Path to parameters file, which includes 
%                spatial bin size and smoothing kernel. 
%                Default = '/Volumes/groups/giocomo/export/data/Projects/ ...
%                JohnKei_NPH3/UniversalParams'

% outputs:
%     all_fr: smoothed firing rate over position for each cell

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
all_fr = nan(nCells, spatialBins); % preallocate matrix of cells' firing rates across spatial bins
fprintf('Calculating firing rate for %d cells with a spatial bin size of %dcm\n',nCells,params.SpatialBin);

%% Calculate firing rates
% specify inputs into the firing rates calculation
trackEnd = trackLength;
p = params;

for k = 1:nCells
    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

    % get spike times and index into post for cell k 
    spike_t = sp.st(sp.clu==cells_to_plot(k));
    
    [~,~,spike_idx] = histcounts(spike_t,post);
    [~,~,spikeTrial_idx] = histcounts(spike_t,trial);
    

    for i = 1:max(trial)
        itrial_kfr = calculateSmoothedFiringRate(spike_idx(i==trial(spike_idx)), posx, p, trackEnd);
        singleCellallTrialsFR(i,:) = itrial_kfr;
%         figure(1)
%         plot(itrial_kfr); %plot FR for single trial in single cell
%         figure(2) % plot rasterplot for single trial in single cell
%         plot(posx(spike_idx(i==trial(spike_idx))),trial(spike_idx(i==trial(spike_idx))),'k.');      
    end
    
%   kfr = calculateSmoothedFiringRate(spike_idx, posx, p, trackEnd); % get average firing rate collapsed across all trials
%   plot(kfr); % plot average firing rate collapsed across all trials

    
    stackedCellVectorTrialsFR = reshape(stackedCellVectorTrialsFR.',1,[]);
    all_fr(k, :) = singleCellallTrialsFR;
    
end

figure();
corrcoefMatrix = corrcoef(singleCellallTrialsFR');
imagesc(corrcoefMatrix);
colorbar;
set(gca,'XTick',0:10:400);
xticklabels(xticks-100)
set(gca,'YTick',0:10:400);
yticklabels(yticks-100)


first = regexp(matPath, 'g0/') + 3;
saveName = strcat(matPath(first: end-4), '_firing rates');

save(saveName, 'all_fr');

% save all_fr in same directory as .mat folder

% parse out the .mat filename from before and use it to save the new .mat
% file.



end