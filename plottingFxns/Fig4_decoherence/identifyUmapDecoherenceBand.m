function [dch_idx,dchTimeSec,dchStartDelaySec] = identifyUmapDecoherenceBand(clusterIdentifiers, cells)
% IDs the decoherence band in number of trials from the results of UMAP being run on a
% spatially-binned FR across all cells

% Find index trial cluster regions
clusterBoundaries = diff(clusterIdentifiers) ~= 0;
clstBegin = find([true,clusterBoundaries])';  % beginning indices
clstEnd = find([clusterBoundaries,true])';    % ending indices
clstLen = 1 + clstEnd - clstBegin;           % cluster length



clst = [clstBegin, clstEnd, clstLen];

% Pull out first cluster of trials after 100 
postKetClst = clst(find(clst(:,1) > 100),:);

% return the indices of the cluster of interest
dch_idx = []; dchTimeSec=[]; dchStartDelaySec=[];
try
   for i = 1:size(postKetClst,1)
       row = postKetClst(i,:);
       dch_idx_temp = row(1):row(2);
       dchTimeSec_temp = findDecoherencePeriodLength(dch_idx_temp, cells);
       dchTimeMin = dchTimeSec_temp/60;

       dchStartDelaySec_temp = findDchPeriodStartDelay(dch_idx_temp, cells);
       dchStartDelayMin = dchStartDelaySec_temp/60;
       
       % return cluster of dchTime is longer than 5 min and if it occurs
       % less than 30 min after the ketamine injection
       if dchTimeMin>5 && dchStartDelayMin<30 
           dch_idx = dch_idx_temp;
           dchTimeSec=dchTimeSec_temp;
           dchStartDelaySec=dchStartDelaySec_temp;
           break
       end
   end
   
catch
   dch_idx = []; dchTimeSec=[]; dchStartDelaySec=[];
end



end

% Keep only cluster of length 10+trials
% clstIdx = find(clstLen > 10);
% clst = [clstBegin(clstIdx), clstEnd(clstIdx), clstLen(clstIdx)];