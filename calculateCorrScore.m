function [drugCorrEffectScore, cellCorrScore, corrTemplate] = calculateCorrScore(singleCellallTrialsFR, trials_corrTemplate, trial)
% calculate correlation score between trial FR compared to baseline template
% and generate a value that represents the effect of a drug on the
% correlation score
% Input: 
%   singleCellallTrialsFR - firing rate maps by trials
%   trials_corrTemplate - the number of trials in baseline template
% Output:
%   drugCorrEffectScore - difference between the peak of the corrScore
%       in the 25 trials before drug and the trough of the corrScore in the
%       25 trials after drug appication (assumes drug is 50 trials after
%       the baseline template)
%   cellCorrScore - correlation score between a trial and the baseline
%       template
%   corrTemplate - generates baseline template from mean FR map of the baseline trials


% for cell k, calculate correlation template from 1st 40trials
corrTemplate = mean(singleCellallTrialsFR(1:trials_corrTemplate,:));
%     plot(corrTemplate)
cellCorrScore = nan(numel(trials_corrTemplate:max(trial)),1);

for i = trials_corrTemplate:max(trial)
    corrMatrix = corrcoef(corrTemplate, singleCellallTrialsFR(i,:));
    corrScore = corrMatrix(1,2);
    cellCorrScore(i-trials_corrTemplate+1) = corrScore;
end

cellCorrScore = fillmissing(cellCorrScore,'spline');

smoothCellCorrScore = smoothdata(cellCorrScore, 'gaussian',10);

corrPeak_preKet = max(smoothCellCorrScore(25:50));
corrTrough_postKet = min(smoothCellCorrScore(50:75));
drugCorrEffectScore = corrPeak_preKet - corrTrough_postKet;

end