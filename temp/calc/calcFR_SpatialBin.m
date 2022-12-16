function [firing_rate, frPerTrial] = calcFR_SpatialBin(idx, trial_posx,posx, p, trackEnd)
% calculates smoothed firing rate on linear track
% Malcolm Campbell 5/21/15
% modified 6/6/18 MGC
% all time bins are now equal length
% edited kei masuda 1/14/20
% inputs:
%     idx: spikes indexs (all spikes with cell IDs attached)
%     posx: position of animal
%     p: params (spatial bin size, etc)
%     TrackEnd: typically, the length of the track 
% outputs:
%     firing_rate: smoothed firing rate over position

% divide spike counts by occupancy


binedges = p.TrackStart:p.SpatialBin:trackEnd;
time_per_bin = histcounts(trial_posx, binedges);
time_per_bin = time_per_bin * p.TimeBin;
spikesPerBin = histcounts(posx(idx), binedges);
firing_rate = spikesPerBin./time_per_bin;
frPerTrial = sum(firing_rate)/sum(time_per_bin);
% % interpolate missing values
if sum(isnan(firing_rate))>0
    %firing_rate(isnan(firing_rate)) = 0;
    try
        firing_rate = interp1(find(~isnan(firing_rate)),firing_rate(~isnan(firing_rate)),1:numel(firing_rate));
    catch
        if sum(~isnan(firing_rate)) < 2 % if there aren't at least 2 data points to be able to interpolate then set FR to 0
            firing_rate(isnan(firing_rate)) = 0;
        end
    end
end


end
