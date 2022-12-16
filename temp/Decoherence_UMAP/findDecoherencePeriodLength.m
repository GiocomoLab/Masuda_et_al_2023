function dchTimeSec = findDecoherencePeriodLength(dch_idx, cells)
% Finds how long in seconds the decoherence period is; assumes the sampling
% rate is 20 Hz
% Input:
% - dch_idx: row vector of trials identified as the decoherence period
% - cells: standard cells struct
% Output:
% - dchTimeSec: double with value of decoherence period in seconds

    
    
    trials = cells.trial(1).trial;
    post = cells.posT(1).post;
    
    Fs = 1/(post(2)-post(1)); %sampling rate is 50hz
    
    trialStart = min(dch_idx);
    trialEnd = max(dch_idx);
    
    try
        idxTrialStart = find(trials==trialStart,1,'first');
        idxTrialEnd = find(trials==trialEnd,1,'last');
        dchTimeSec = (idxTrialEnd-idxTrialStart)/Fs;
    catch
        dchTimeSec = 0;
    end
    
end