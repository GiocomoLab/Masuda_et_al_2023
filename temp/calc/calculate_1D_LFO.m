function linearFractionalOccupancy = calculate_1D_LFO(posx,post,spatialBin, speed)

    trackLength = 400;
    trackStart = 0;
    posSamplingRate = mean(diff(post)); % sec per frame
    
    speedThreshold = 2; %2cm/s
    [~,spdFltrPosX,spdFltrPosT] = speedFilterData([],posx, post, speed,speedThreshold);
    
    binedges = trackStart:spatialBin:trackLength;
    time_per_bin = histcounts(spdFltrPosX, binedges);
    time_per_bin = time_per_bin * posSamplingRate;
    total_time = numel(spdFltrPosT) * posSamplingRate;
    linearFractionalOccupancy = time_per_bin/total_time;
    
end