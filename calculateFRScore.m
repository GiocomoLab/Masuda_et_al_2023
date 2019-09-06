function drugFREffectScore = calculateFRScore(singleCellallTrialsFR, preDrugTrials)
% calculate a score that represents the effect of drug on FR
% Input: 
%   singleCellallTrialsFR - firing rate maps by trials
%   preDrugTrials - the number of trials before drug is administered
% Output:
%   drugFREffectScore - difference between the avg FR of the preDrug Trials
%       and the the avg FR of the equivalent number of trials post drug
%       administration

trialAvgFR = mean(singleCellallTrialsFR,2);
preDrugFR = mean(trialAvgFR(1:preDrugTrials));
postDrugFR = mean(trialAvgFR(preDrugTrials:preDrugTrials*2));
drugFREffectScore = postDrugFR - preDrugFR;

end