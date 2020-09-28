function pool_decoherence_sessions = poolDecoherenceSessions(cells)
% pool decoherence sessions from cells struct
% can be used to pool together decoherence 

seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));
pool_decoherence_sessions(numel(seshes)) = struct(); 

for i = 1:numel(seshes)
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    cellIndx = find(seshIndx, 1, 'first');
    
    pool_decoherence_sessions(i).name = cells.metadata{cellIndx,2};
    pool_decoherence_sessions(i).genotype = cells.metadata{cellIndx,4};
    pool_decoherence_sessions(i).sessionDate = cells.metadata{cellIndx,3};
    
    
    pool_decoherence_sessions(i).reduction = cells.dch(cellIndx).dch.umapOutput.reduction;
    pool_decoherence_sessions(i).clusterIdentifiers = cells.dch(cellIndx).dch.umapOutput.clusterIdentifiers;
    
    pool_decoherence_sessions(i).idxCellArray = cells.dch(cellIndx).dch.idxCellArray;
    pool_decoherence_sessions(i).idxClusterArray = cells.dch(cellIndx).dch.idxClusterArray;
    pool_decoherence_sessions(i).decoherenceIdx = cells.dch(cellIndx).dch.decoherenceIdx;
    pool_decoherence_sessions(i).seshName = cells.dch(cellIndx).dch.seshName;
    pool_decoherence_sessions(i).decoherenceTime = cells.dch(cellIndx).dch.decoherenceTime;
    pool_decoherence_sessions(i).decoherenceStartDelay = cells.dch(cellIndx).dch.decoherenceStartDelay;
    pool_decoherence_sessions(i).decoherenceTimeIdx = cells.dch(cellIndx).dch.decoherenceTimeIdx;
    pool_decoherence_sessions(i).Fs = cells.dch(cellIndx).dch.Fs;
    pool_decoherence_sessions(i).umapOutput = cells.dch(cellIndx).dch.umapOutput;
    % uncomment this next time I run calculate decoherence sessions
%     pool_decoherence_sessions(i).timeDownSample = cells.dch(cellIndx).dch.timeDownSample; 
    pool_decoherence_sessions(i).timeDownSample = 100;
end

end