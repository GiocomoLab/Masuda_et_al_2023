function nnMatrix(matPath)

addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));
% 
% load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/G4_190620_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');
% load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/HCN1_190619_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');
% load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/G3_190708_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');

load(matPath);
%%
cells_to_plot = sp.cids(sp.cgs==2); 
nCells = size(cells_to_plot, 2);

% compute some useful information (like spike depths)
[~, spikeDepths, ~, ~, ~, ~, ~] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

% get spike depths
spike_depth = nan(numel(cells_to_plot),1);
for k = 1:numel(cells_to_plot)
    spike_depth(k) = median(spikeDepths(sp.clu==cells_to_plot(k)));
end

% sort cells_to_plot by spike_depth (descending)
[spike_depth,sort_idx] = sort(spike_depth,'descend');
cells_to_plot = cells_to_plot(sort_idx);

%% Calculating firing rate

% calculate firing rate by time
all_fr = []; 
ds_factor = 50;
smoothSigma = 10;
trial_ds = downsample(trial, ds_factor); 
for i = 1:nCells
   
    % get spike times for cell i
    spike_t = sp.st(sp.clu==cells_to_plot(i));
    
    % calculate firing rate
    fr = calcSmoothedFR_Time(post, spike_t, ds_factor, smoothSigma);
%     fr = calcSmoothedFR_SpatialBin(idx, posx, p, TrackEnd);
    all_fr(i, :) = fr;
        
end

all_fr = all_fr';

%%
nn_matrix = corr(all_fr(trial_ds>0 & trial_ds<=300,:));

test_matrix = corr(all_fr(trial_ds>110 & trial_ds<=150,:));
% sort the rows/columns
% generate distances between each row
D = pdist(nn_matrix);

% create tree
tree = linkage(D,'average');

% find the order
leafOrder = optimalleaforder(tree,D,'transformation','inverse');

% re-sort the matrix
sorted_matrix = nn_matrix(leafOrder,:);
sorted_matrix = sorted_matrix(:,leafOrder);

% re-sort the test matrix
sorted_testMatrix = test_matrix(leafOrder,:);
sorted_testMatrix = sorted_testMatrix(:,leafOrder);

% plot
figure(1)
subplot(1,3,1)
imagesc(nna_matrix)
title('original baseline matrix')

subplot(1,3,2)
imagesc(sorted_matrix);
title('clustered baseline matrix matrix')

subplot(1,3,3)
imagesc(sorted_testMatrix);
title('sorted test matrix')
end