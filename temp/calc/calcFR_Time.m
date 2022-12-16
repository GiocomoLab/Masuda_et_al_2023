function firing_rate = calcFR_Time(post, spike_t, ds_factor)

% Calculate firing rate from spike times for one cell
% INPUTS
%   post (n x 1): sampling times. E.g., [0; 0.02; 0.04 ...]
%   spike_t (m x 1): the l spike times for a cell. E.g., [1.13; 1.16; 1.18]
%   ds_factor (scalar): indicates how much to downsample post. E.g., 20.
%
% OUTPUTS
%   firing_rate (l x 1): firing rate over time. Number of time bins is a
%       function of how much downsampling is used (l = ceil(m/ds_factor)).
%
% Created by John Wen 190712
% Last edited by FKM 8/28/19

%% Downsample post
post_ds = downsample(post, ds_factor); % downsample to count spikes in timebin
timebin = diff(post_ds(1:2)); 

% append timebin to catch any spike times that fall outside range after downsampling
post_ds = [post_ds; max(post_ds) + timebin]; 
%% Calculate timebin
% timebin = diff(post(1:2)); 
% 
% % append timebin to catch any spike times that fall outside range after downsampling
% post = [post; max(post) + timebin]; 

%% Bin spike times into post and calculate smoothed firing rate
firing_rate = histcounts(spike_t, post_ds)/timebin; % divide by time to get fr


%% Downsample firing rate
% firing_rate = downsample(firing_rate, ds_factor);

end