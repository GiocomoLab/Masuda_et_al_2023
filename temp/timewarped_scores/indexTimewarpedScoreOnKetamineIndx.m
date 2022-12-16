function ketIndx_timewarpedScore = indexTimewarpedScoreOnKetamineIndx(tw,min)
% take timewarped scores and chop everything before ketamine off

if ~exist('paramsPath','var')
    params = readtable('./UniversalParams.xlsx');
else
    params = readtable(paramsPath);
end

samplingRate = params.TimeBin;
ds_factor = params.ds_factor;
secsInMin = 60; 
conversionFactor = secsInMin/(ds_factor*samplingRate);
conversionRate = min * conversionFactor;

A = tw.timewarpedScore;
ketIndx = tw.ketamineIndx;

if ~isempty(min)
    endIndx = ketIndx+conversionRate;
    for k=1:numel(A)
      A{k}=A{k}(ketIndx(k):endIndx(k));
    end 
elseif ~isnan(tw.gainIndx)
    endIndx = tw.gainIndx;
    for k=1:numel(A)
      A{k}=A{k}(ketIndx(k):endIndx(k));
    end
else
    for k=1:numel(A)
      A{k}=A{k}(ketIndx(k):end);
    end
end
    
ketIndx_timewarpedScore = A;