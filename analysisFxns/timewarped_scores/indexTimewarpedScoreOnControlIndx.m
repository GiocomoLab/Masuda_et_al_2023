function cntrlIndx_timewarpedScore = indexTimewarpedScoreOnControlIndx(tw,min)
% take timewarped scores and chop everything before cntrl and after ketamine off

if ~exist('paramsPath','var')
    params = readtable('./UniversalParams.xlsx');
else
    params = readtable(paramsPath);
end

A = tw.timewarpedScore;
ketIndx = tw.ketamineIndx;
cntrlIndx = tw.controlIndx;

samplingRate = params.TimeBin;
ds_factor = params.ds_factor;
secsInMin = 60; 
conversionFactor = secsInMin/(ds_factor*samplingRate);
conversionRate = min * conversionFactor;

if ~isempty(min)
    endIndx = cntrlIndx+conversionRate;
    for k=1:numel(A)
      A{k}=A{k}(cntrlIndx(k):endIndx(k));
    end 
else
    for k=1:numel(A)
      A{k}=A{k}(cntrlIndx(k):ketIndx(k));
    end
end
    
cntrlIndx_timewarpedScore = A;