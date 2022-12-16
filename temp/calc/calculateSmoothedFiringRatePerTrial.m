function firing_rate = calculateSmoothedFiringRatePerTrial(idx, posx, post, trial, p, TrackEnd)
% calculates smoothed firing rate on linear track for all trials for one 
% cell 

% Modified from Malcolm Campbell's calculateSmoothedFiringRate.m
% John Wen 190702

% Calculate firing rate by dividing spike counts by occupancy. 
% inputs:
%     idx: spike times for one cell indexed into post
%     posx: position along track at every timestamp 
%     post: time value since starting Unity session
%     trial: trial number at every timestamp
%     p: params (spatial bin size, etc)
%     TrackEnd: typically, the length of the track 
% outputs:
%     firing_rate: smoothed firing rate over position for each trial



% get some parameters and preallocate matrices
maxTrial = max(trial);

frPerTrial = nan(maxTrial, TrackEnd/p.SpatialBin);

% create gaussian filter for smoothing later
smoothSigma = p.SmoothSigmaFR/p.SpatialBin;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma)

binedges = p.TrackStart:p.SpatialBin:TrackEnd; % define spatial bin edges

for h = 1:maxTrial

    % index for demarcating trials (last trial index)
    ltri = 0;
    
    % calculate time spent per spatial bin on trial h
    posx_h = posx(trial == h); % get the posx for trial h
    time_per_bin = histcounts(posx_h, binedges); % bin posx_h by binedges
    time_per_bin = time_per_bin * p.TimeBin; % multiply time_per_bin by TimeBin to get actual times
    
    % calculate number of spikes per spatial bin on trial h
    idx_h = idx(idx > ltri && idx <= ltri + size(posx_h, 1));
    ltri = ltri + size(posx_h, 1); 
    
    
    idx_h = idx(trial == h); % get all the spikes for this cell on trial h
    firing_rate = histcounts(posx_h(idx_h), binedges); % count spikes per spatial bin on this trial
    firing_rate = firing_rate./time_per_bin; % calculate firing rate for this cell by dividing by occupancy
    
    % interpolate missing values
    if sum(isnan(firing_rate))>0
        firing_rate = interp1(find(~isnan(firing_rate)),firing_rate(~isnan(firing_rate)),1:numel(firing_rate));
    end

    % smooth firing rate
    firing_rate_smooth = conv(repmat(firing_rate,1,3),gauss_filter,'same');
    firing_rate = firing_rate_smooth(numel(firing_rate)+1:numel(firing_rate)*2);
    
    frPerTrial(h, :) = firing_rate; % save firing rate for trial h in to frPerTrial
    
end


end
