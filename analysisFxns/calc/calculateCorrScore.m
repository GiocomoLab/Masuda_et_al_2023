function [drugCorrEffectScore, cellCorrScore, corrTemplate] = calculateCorrScore(singleCellallTrialsFR, trials_corrTemplate)
% calculate correlation score between trial FR compared to baseline template
% and generate a value that represents the effect of a drug on the
% correlation score
% Input: 
%   singleCellallTrialsFR - firing rate maps by trials
%   trials_corrTemplate - the number of trials in baseline template
% Output:
%   drugCorrEffectScore - difference between the average of the corrScore
%       in the 25 trials before drug and the average of the corrScore in the
%       25 trials after drug appication (assumes drug administered 50 trials after
%       the baseline template)
%   cellCorrScore - correlation score between a trial and the baseline
%       template
%   corrTemplate - generates baseline template from mean FR map of the baseline trials


% for cell k, calculate correlation template from 1st 50 trials
corrTemplate = mean(singleCellallTrialsFR(1:trials_corrTemplate,:));
%     plot(corrTemplate)
numTrial = size(singleCellallTrialsFR,1);
% cellCorrScore = nan(numel(trials_corrTemplate:numTrial),1);
cellCorrScore = nan(numel(1:numTrial),1);

% for i = trials_corrTemplate+1:numTrial
for i = 1:numTrial
    corrMatrix = corrcoef(corrTemplate, singleCellallTrialsFR(i,:));
    corrScore = corrMatrix(1,2);
    cellCorrScore(i) = corrScore;
%     cellCorrScore(i-trials_corrTemplate) = corrScore;
end


% smoothCellCorrScore = smoothdata(cellCorrScore, 'gaussian',10);

corrPeak_preKet = nanmean(cellCorrScore(76:100));
corrTrough_postKet = nanmean(cellCorrScore(101:125));
drugCorrEffectScore = corrPeak_preKet - corrTrough_postKet;

end