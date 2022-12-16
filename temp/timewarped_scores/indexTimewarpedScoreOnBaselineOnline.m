function baselineIndx_timewarpedScore = indexTimewarpedScoreOnBaselineOnline(tw,min)
% take timewarped scores and chop everything before ketamine off

if ~exist('paramsPath','var')
    params = readtable('./UniversalParams.xlsx');
else
    params = readtable(paramsPath);
end

A = tw.timewarpedScore;
cntrlIndx = tw.controlIndx;

samplingRate = params.TimeBin;
ds_factor = params.ds_factor;
secsInMin = 60; 
conversionFactor = secsInMin/(ds_factor*samplingRate);
conversionRate = min * conversionFactor;

if ~isempty(min)
    endIndx = 1+conversionRate;
    for k=1:numel(A)
      A{k}=A{k}(1:endIndx);
    end 
else
    for k=1:numel(A)
      A{k}=A{k}(1:cntrlIndx(k));
    end
end

baselineIndx_timewarpedScore = A;