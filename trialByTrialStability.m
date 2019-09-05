function rho = trialByTrialStability(testTrials, trial, spike_idx, posx, p, trackEnd)
% calculate stability comparing every two trials and then averaging rho
% Inputs
% -testTrials: range of trials e.g. 1:50
%    * if range is odd, then last trial is correlated to previous trial
% -trial: vector with trial info by sample
% -spike_idx: spikes indexed into unity samples
% -posx: vector of unity positions
% -p: universal parameters
% -trackEnd: length of VR track
% Outputs
%- rho: stability value
    
    numAvgTrials = 5;
    allTestRho = nan(numel(min(testTrials):numAvgTrials:max(testTrials)),1);
    for i = min(testTrials):numAvgTrials*2:max(testTrials)
        fr1 = calcSmoothedFR_SpatialBin(spike_idx(trial(spike_idx)>=i & trial(spike_idx)<(i+numAvgTrials)), posx(trial>=i & trial<(i+numAvgTrials)),posx, p, trackEnd);
        j = i+numAvgTrials;
        if j > max(testTrials)
            fr2 = calcSmoothedFR_SpatialBin(spike_idx(trial(spike_idx)>=j & trial(spike_idx)<(j-numAvgTrials)), posx(trial>=j & trial<(j-numAvgTrials)),posx, p, trackEnd);
        else
            fr2 = calcSmoothedFR_SpatialBin(spike_idx(trial(spike_idx)==i+1), posx(trial==i+1),posx, p, trackEnd);
        end

        testRho = corr(fr1',fr2');
        allTestRho(floor(i/2)+1) = testRho;
    end
    rho = mean(allTestRho);
    
    
end