function [cells_to_plot, spike_depth,waveforms] = spPreProcess(sp)
% Input: 
%     sp - struct
% Output: 
%   cells_to_plot - cluster IDs of good cells sorted by spikeDepth
%   from most dorsal to most ventral on the probe
%   spike_depth - median spike depth for each cell
%   waveforms - template waveforms for all the clusters 

addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
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
    
end