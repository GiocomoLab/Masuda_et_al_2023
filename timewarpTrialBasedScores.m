function tw = timewarpTrialBasedScores(cells, trialBasedScore,ds_factor)

trial = cells.trial(1).trial;
ds_trial = downsample(trial,ds_factor);

timewarpedScore = nan(size(trialBasedScore,1),size(ds_trial,1));

    for i = 1:size(timewarpedScore,2)
       timewarpedScore(:,i) = trialBasedScore(:,ds_trial(i));
    end

tw.timewarpedScore = timewarpedScore;

ketamineTrial = 100;
controlTrial = 50;
tw.ketamineIndx = find(ds_trial > ketamineTrial,1,'first');
tw.controlIndx = find(ds_trial > controlTrial,1,'first');
    
end
