function ketIndx_timewarpedScore = indexTimewarpedScoreOnKetamineIndx(tw)
% take timewarped scores and chop everything before ketamine off

A = tw.timewarpedScore;
ketIndx = tw.ketamineIndx;
gainIndx = tw.gainIndx;

for k=1:numel(A)
  A{k}=A{k}(ketIndx(k):gainIndx(k));
end

ketIndx_timewarpedScore = A;