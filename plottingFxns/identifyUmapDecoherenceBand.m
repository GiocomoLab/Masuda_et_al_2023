function dch_idx = identifyUmapDecoherenceBand(clusterIdentifiers)

% Find index trial cluster regions
clusterBoundaries = diff(clusterIdentifiers) ~= 0;
clstBegin = find([true,clusterBoundaries])';  % beginning indices
clstEnd = find([clusterBoundaries,true])';    % ending indices
clstLen = 1 + clstEnd - clstBegin;           % cluster length

% Keep only cluster of length 10+
clstIdx = find(clstLen > 10);
clst = [clstBegin(clstIdx), clstEnd(clstIdx), clstLen(clstIdx)];

% Pull out first cluster of trials after 100 that's greater than 5 trials
row = clst(find(clst(:,1) > 100, 1),:);
% val = clusterIdentifiers(row(1));

% return the indices of the cluste rof interest
try
   dch_idx = row(1):row(2);
catch
   dch_idx = []; 
end
end