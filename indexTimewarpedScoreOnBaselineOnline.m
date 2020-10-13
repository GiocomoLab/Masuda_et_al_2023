function baselineIndx_timewarpedScore = indexTimewarpedScoreOnBaselineOnline(tw)
% take timewarped scores and chop everything before ketamine off

A = tw.timewarpedScore;
cntrlIndx = tw.controlIndx;

for k=1:numel(A)
  A{k}=A{k}(1:cntrlIndx(k));
end

baselineIndx_timewarpedScore = A;