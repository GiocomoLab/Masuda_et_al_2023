function stabilityScore = calculateStabilityScore(singleCellallTrialsFR)
% first normalizes the firing rate by centering the data to have mean 0
% calculate a stability score between trial FR compared to the trial before
% and the trial after it. If it's the first trial compare with just trial
% after it, and if it's the last trial compare with just the trial before
% it
% 
% Input: 
%   singleCellallTrialsFR - firing rate maps x trials
% Output:
%   stabilityScore - singleCellallTrialsFR
    
    zscore_FR = normalize(singleCellallTrialsFR,2,'center');
    imagesc(zscore_FR)
    stabilityScore = nan(size(singleCellallTrialsFR,1),1);
    for i = 1:numel(stabilityScore)
        
        if i == 1
            fr1 = singleCellallTrialsFR(i,:); 
            fr2 = singleCellallTrialsFR(i+1,:);
            rho = corr(fr1',fr2');
            
        elseif i == numel(stabilityScore)
            fr1 = singleCellallTrialsFR(i,:); 
            fr2 = singleCellallTrialsFR(i-1,:);
            rho = corr(fr1',fr2');
        else
            fr0 = singleCellallTrialsFR(i-1,:); 
            fr1 = singleCellallTrialsFR(i,:); 
            fr2 = singleCellallTrialsFR(i+1,:);
            rho1 = corr(fr1',fr2');
            rho2 = corr(fr0',fr1');
            rho = mean([rho1,rho2]);
        end
      
        stabilityScore(i) = rho;
    end
    

end