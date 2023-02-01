my_cells = readtable('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/F1/F1_190620_johncontrasttrack_train1_sp_cells_by_eye.xlsx');
my_cells = table2array(my_cells);

sp_cells_st = struct('cell_number', [], 'fr', [], 'MI', [], 'early_trial_normalized_correlation', []);

sp_cells_st.cell_number = my_cells;


for i = 1:size(my_cells, 1)
    
    % use cell number to find in cells_to_plot the corresponding MI value
    % and cell depth
    cells_to_plot_idx = find(cells_to_plot == my_cells(i));
    sp_cells_st.MI = [sp_cells_st.MI MI_vec(cells_to_plot_idx)];
    
    % use cell number to calculate spatial firing rate and calculate
    % correlation then put in struct too
    
end
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


%% Preallocate and store data for all cell's firing rates
nCells = size(my_cells, 1);
spatialBins = trackLength/2; 
all_fr = nan(nCells, max(trial),spatialBins); % preallocate matrix of cells' firing rates across spatial bins

all_corrmatrix = nan(nCells, max(trial), max(trial)); % preallocate matrix of cells' correlations between trials
fprintf('Calculating firing rate for %d cells with a spatial bin size of %dcm\n',nCells,params.SpatialBin);

%% Calculate firing rates
% specify inputs into the firing rates calculation
trackEnd = trackLength;
p = params;

% calculate the firing rate for a single cell across all trials

for k = 1:nCells
    fprintf('cell %d (%d/%d)\n',my_cells(k),k,numel(my_cells));

    % get spike times and index into post for cell k 
    spike_t = sp.st(sp.clu==my_cells(k));

    [~,~,spike_idx] = histcounts(spike_t,post);
    [~,~,spikeTrial_idx] = histcounts(spike_t,trial); % NOT SURE IF WE STILL NEED THIS


    % for cell k, iteratively calculate the firing rate for each trial
    for i = 1:max(trial)
        itrial_kfr = calculateSmoothedFiringRate(spike_idx(i==trial(spike_idx)), posx, p, trackEnd);
        singleCellallTrialsFR(i,:) = itrial_kfr;
%         figure(1)
%         plot(itrial_kfr); %plot FR for single trial in single cell
%         figure(2) % plot rasterplot for single trial in single cell
%         plot(posx(spike_idx(i==trial(spike_idx))),trial(spike_idx(i==trial(spike_idx))),'k.');      
    end

    % calculate trial by trial correlations for one cell 
    corrMatrix = corr(singleCellallTrialsFR');


    % store one cell's firing rates and trial by trial correlation matrix
    % in a session matrix containing all cells
    all_fr(k, :, :) = singleCellallTrialsFR;
    all_corrmatrix(k, :, :) = corrMatrix;

end

sp_cells_st.fr = all_fr;

%% Calculate Pearson correlation for top 10 trials compared to average

avg_fr_mat = squeeze(mean(all_fr, 2)); % gives trial-averaged fr across spatial bins

for i = 1:nCells
    
    % average first 10 trials of firing
    early = squeeze(mean(all_fr(i, 1:10, :), 2));
    
    % average trials 50-100
    mid = squeeze(mean(all_fr(i, 50:100, :), 2));
    
    % calculate Pearson correlations for beginning and middle trials
    early_corr = corr(early, avg_fr_mat(i, :)');
    mid_corr = corr(mid, avg_fr_mat(i, :)');
    
    % take ratio and store in struct
    early_corr_normalized = early_corr/mid_corr;
    sp_cells_st.early_trial_normalized_correlation = ...
        [sp_cells_st.early_trial_normalized_correlation early_corr_normalized];
    
end

tbl = table(sp_cells_st.early_trial_normalized_correlation', sp_cells_st.MI');
tbl.Properties.VariableNames = {'Early_Trial_Correlations', 'MI'};
mdl = fitlm(tbl, 'linear');



