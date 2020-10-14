function tw = timewarpTrialBasedScores(cells, trialBasedScore)


if ~exist('paramsPath','var')
    params = readtable('./UniversalParams.xlsx');
else
    params = readtable(paramsPath);
end

ds_factor = params.ds_factor;
ketamineTrial = params.ketamineTrial;
controlTrial = params.controlTrial;
gainTrial = params.gainTrial;

trial = cells.trial;

numCells = size(trial,2);
ketamineIndx = nan(numCells,1);
controlIndx = nan(numCells,1);
gainIndx = nan(numCells,1);
timewarpedScores = cell(numCells,1);

for i = 1:numCells
    
    ds_trial = downsample(trial(i).trial,ds_factor);
    numSamples = size(ds_trial,1);
    
    ketamineIndx(i) = find(ds_trial > ketamineTrial,1,'first');
    controlIndx(i) = find(ds_trial > controlTrial,1,'first');
    
    if max(ds_trial) > gainTrial
        gainIndx(i) = find(ds_trial > gainTrial,1,'first');
    else
        gainIndx(i) = nan;
    end
    
    
    single_cell_trialBasedScore = trialBasedScore(i,:);
    single_cell_timewarpedScore = nan(numSamples,1);
    for j = 1:numSamples
        try
            single_cell_timewarpedScore(j) = single_cell_trialBasedScore(ds_trial(j));
        catch
            if ds_trial(j) > 300
            	single_cell_timewarpedScore(j) = single_cell_trialBasedScore(300);
            else
                
            end
        end
    end
    
    timewarpedScores{i} = single_cell_timewarpedScore;
end

tw.timewarpedScore = timewarpedScores;
tw.ketamineIndx = ketamineIndx;
tw.controlIndx = controlIndx;
tw.gainIndx = gainIndx;

end
