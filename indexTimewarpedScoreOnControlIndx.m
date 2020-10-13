function cntrlIndx_timewarpedScore = indexTimewarpedScoreOnControlIndx(tw)
% take timewarped scores and chop everything before cntrl and after ketamine off

A = tw.timewarpedScore;
ketIndx = tw.ketamineIndx;
cntrlIndx = tw.controlIndx;

for k=1:numel(A)
  A{k}=A{k}(cntrlIndx(k):ketIndx(k));
end

cntrlIndx_timewarpedScore = A;