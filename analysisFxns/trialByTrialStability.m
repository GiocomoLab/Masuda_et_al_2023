function rho = trialByTrialStability(frMap, testTrials, numAvgTrials)
% calculate stability comparing every two trials and then averaging rho
% Inputs
% - frMap - trial x spatialBin firing rate map
% - testTrials: range of trials e.g. 1:50
%    * if range is odd, then last trial is correlated to previous trial
% - numAvgTrials - test stability between first and second half of every N contiguous trials must be even e.g. 10
% Outputs
%- rho: stability value
    
    if ~mod(numAvgTrials, 2)
        allTestRho = nan(numel(min(testTrials):numAvgTrials:max(testTrials)),1);
        for i = min(testTrials):numAvgTrials:max(testTrials)
            numHalf = numAvgTrials/2;
            fr1 = frMap(i:i+(numHalf-1),:); %get first half 
            fr1 = reshape(fr1.',1,[]); % flatten so trials are concatenated
            
            fr2 = frMap(i+numHalf:i+(2*numHalf-1),:);
            fr2 = reshape(fr2.',1,[]); % flatten so trials are concatenated
            
            testRho = corr(fr1',fr2');
            index = ((i-min(testTrials))/numAvgTrials)+1;
            allTestRho(index) = testRho;
        end
        rho = nanmean(allTestRho);
    else
       fprintf('\nnumAvgTrials must be even');
    end
    
end