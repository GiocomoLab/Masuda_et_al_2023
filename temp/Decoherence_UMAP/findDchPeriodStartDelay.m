function dchStartDelaySec = findDchPeriodStartDelay(dch_idx, cells)
% Finds how long in seconds the decoherence period is; assumes the sampling
% rate is 20 Hz
% Input:
% - dch_idx: row vector of trials identified as the decoherence period
% - cells: standard cells struct
% Output:
% - dchStartDelaySec: double with value of decoherence period in seconds

    
    
    trials = cells.trial(1).trial;
    post = cells.posT(1).post;
    
    Fs = 1/(post(2)-post(1)); %sampling rate is 50hz
    
    trialStart = min(dch_idx);
    ketamineStart = 101;
    try
        idxTrialStart = find(trials==trialStart,1,'first');
        idxKetamineStart = find(trials==ketamineStart,1,'first');
        dchStartDelaySec = (idxTrialStart-idxKetamineStart)/Fs;
    catch
        dchStartDelaySec = 0;
    end
    
end