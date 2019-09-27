function [drugFRdiff,cntrlFRdiff,drugFREffectScore] = calculateFRScore(singleCellallTrialsFR, preDrugTrials, preControlTrials)
% calculate a score that represents the effect of drug on FR
% Input: 
%   singleCellallTrialsFR - firing rate maps by trials
%   preDrugTrials - the number of trials before drug is administered
%   preControlTrials - the number of trials before control injx is administered
% Output:
%   drugFREffectScore - the ratio drug vs cntrl of the differences between the avg FR of the preInjx Trials
%       and the the avg FR of the equivalent number of trials postInjx

trialAvgFR = nanmean(singleCellallTrialsFR,2);
preDrugFR = nanmean(trialAvgFR(1:preDrugTrials));
postDrugFR = nanmean(trialAvgFR(preDrugTrials:preDrugTrials*2));

preCntrlFR = nanmean(trialAvgFR(1:preControlTrials));
postCntrlFR = nanmean(trialAvgFR(preControlTrials:preDrugTrials));

drugFRdiff = postDrugFR-preDrugFR;
cntrlFRdiff = postCntrlFR-preCntrlFR;

drugFREffectScore = drugFRdiff/cntrlFRdiff;

end