
function dchTimeSec = findDecoherencePeriodLength(dch_idx, cells)
    Fs = 0.05;
    
    trials = cells.trial(1).trial;
    post = cells.posT(1).post;
    
    trialStart = min(dch_idx);
    trialEnd = max(dch_idx);
    idxTrialStart = find(trials==trialStart,1,'first');
    idxTrialEnd = find(trials==trialEnd,1,'last');
    dchTimeSec = (idxTrialEnd-idxTrialStart)*Fs;
    
end