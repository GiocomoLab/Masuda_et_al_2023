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

% ds_trial = downsample(trial(i).trial,ds_factor);
% timeIndexSize = size(ds_trial,1);
% timewarpedScore = nan(size(trialBasedScore,1),timeIndexSize);
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
    gainIndx(i) = find(ds_trial > gainTrial,1,'first');

    single_cell_trialBasedScore = trialBasedScore(i,:);

    % ketamineAlignedTimeIndexSize = timeIndexSize - tw.ketamineIndx;
    % aligned_to_ketamine_start_timewarpedScore = nan(size(trialBasedScore,1),ketamineAlignedTimeIndexSize); 
    single_cell_timewarpedScore = nan(numSamples,1);
    for j = 1:numSamples
        single_cell_timewarpedScore(j) = single_cell_trialBasedScore(ds_trial(j));
    end
%    if ds_trial(i) > 100
%        aligned_to_ketamine_start_timewarpedScore
%    end
    
    timewarpedScores{i} = single_cell_timewarpedScore;
end

tw.timewarpedScore = timewarpedScores;

end
