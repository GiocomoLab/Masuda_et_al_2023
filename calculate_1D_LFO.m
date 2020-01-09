function linearFractionalOccupancy = calculate_1D_LFO(posx,post)
    if ~exist('paramsPath','var')
        addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
        params = readtable('UniversalParams.xlsx');
    else
        params = readtable(paramsPath);
    end

    trackLength = floor(max(posx));
    trackStart = 0;
    posSamplingRate = mean(diff(post)); % sec per frame
    

    binedges = trackStart:params.SpatialBin:trackLength;
    time_per_bin = histcounts(posx, binedges);
    time_per_bin = time_per_bin * posSamplingRate;
    total_time = max(post);
    linearFractionalOccupancy = time_per_bin/total_time;
end