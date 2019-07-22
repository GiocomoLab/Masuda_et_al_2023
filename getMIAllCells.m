function [MI_vec, MI_th, all_fr, spatial_cells] = ...
    getMIAllCells(matPath, ds_factor, smoothSigma, fr_bins, spatial_binsize, trackLength)

% Calculate mutual information for all cells in file specified by matPath
%
% INPUTS
%   matPath: path to .mat file after running sync_vr_to_np.m. Specify as
%       string. e.g.
%       '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/E1/...
%       E1_190617_johncontrasttrack_train1_g0/E1_190617_johnrepeatingtrack_train1.mat'
%   ds_factor (scalar): indicates how much to downsample post to calculate 
%       firing rate. E.g., 10.
%   smoothSigma (scalar): indicate degree of firing rate smoothing. E.g., 10.
%   fr_bins (scalar): indicate number of firing rate bins for mutual
%       information calculation. E.g., 20.
%   spatial_binsize (scalar): indicate size of spatial bins for mutual
%       information calculation. E.g., 2.
%   trackLength (scalar): length of track. E.g., 320.
%
% OUTPUTS
%   MI_vec (cells x 1): MI values for each cell 
%   MI_th (cells x 1): MI thresholds for each cell (95% on shuffled data)
%   all_fr (cells x ceil(size(post, 1)/ds_factor)): fr of each cell across
%       time
%   spatial_cells (m x 1): cell numbers of spatial cells
%   
%
% Created by John Wen 190712
% Last edited by John Wen 190712

%% Set number of shuffles to calculate MI significance
nShuffles = 100;

%% Load data, add paths, sort cells by depth
% load specific components from .mat file
load(fullfile(matPath), 'post','posx','sp','trial'); 

addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));

% calculate MI for only good cells
cells_to_plot = sp.cids(sp.cgs==2); 
nCells = size(cells_to_plot, 2);

% compute some useful information (like spike depths)
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

% get spike depths
spike_depth = nan(numel(cells_to_plot),1);
for k = 1:numel(cells_to_plot)
    spike_depth(k) = median(spikeDepths(sp.clu==cells_to_plot(k)));
end

% sort cells_to_plot by spike_depth (descending)
[spike_depth,sort_idx] = sort(spike_depth,'descend');
cells_to_plot = cells_to_plot(sort_idx);

%% Calculate MI for each cell
MI_vec = nan(nCells, 1);
% all_fr = nan(nCells, ceil(size(posx, 1)/ds_factor)); %CHANGE BACK
all_fr = [];
MI_th = nan(nCells, 1);

for i = 1:nCells
    % get spike_t for cell i 
    spike_t = sp.st(sp.clu==cells_to_plot(i));
    
%     % optional code to keep only the spikes up until a certain trial number
%     [~,~,spike_idx] = histcounts(spike_t,post); % index spike_t into post
%     spike_t = spike_t(140 >= trial(spike_idx)); % keep only the spike times happening within 140 trials
%     post_sub = post(trial <= 140); % create new post within 140 trials
%     posx_sub = posx(trial <= 140); % create new posx within 140 trials

    % calculate firing rate
    fr = calculateSmoothFRbyTime(post, spike_t, ds_factor, smoothSigma); 
    all_fr(i, :) = fr;
    max_fr = max(fr);
    
    % downsample posx to match firing rate
    posx_ds = downsample(posx', ds_factor);
    
    % create X_edges and Y_edges for MI calculation
    fr_binsize = max_fr/(fr_bins - 1);
    X_edges = linspace(0, max_fr + fr_binsize, fr_bins);
    
    Y_edges = [0:spatial_binsize:trackLength];
    
    % calculate MI for cell i
    [MI_cell, ~, ~] = getMI(fr, posx_ds, X_edges, Y_edges);
    MI_vec(i) = MI_cell;
    
   
    % determine if MI for cell i is significant (not saving the shuffled
    % distributions though)
    MI_shuffled = nan(nShuffles, 1);
    
    for j = 1:nShuffles
        
        phase = randi(size(posx_ds, 2), 1);
        shuffle_fr = circshift(fr, phase); % circularly shuffle data
        
        % calculate MI for shuffled fr  
        MI_shuffled(j) = getMI(shuffle_fr, posx_ds, X_edges, Y_edges);   
        
    end
    
    MI_th(i) = prctile(MI_shuffled, 95);
    
end

% get the cells that have MI values greater than 95 percentile and are in
% the top 50 percentile for MI 
spatial_cells = cells_to_plot(MI_vec > MI_th & MI_vec > prctile(MI_vec, 50)); 



end




