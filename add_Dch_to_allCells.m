function add_Dch_to_allCells(allCells,dch)
% Find unique sessions in cells struct
seshes = unique(cellfun(@num2str,allCells.metadata(:,1),'uni',0));
count = numel(allCells.metadata(:,1));

clear allCells_Dch
allCells_Dch(count) = struct(); 
z = 0;
fprintf('Done with the Pre-Allocation\n');

for i = 1:numel(seshes)
    clear dchSeshOnly
    dchSeshOnly.idxCellArray = dch.idxCellArray{i};
    dchSeshOnly.idxClusterArray = dch.idxClusterArray{i};
    dchSeshOnly.decoherenceIdx = dch.decoherenceIdx{i};
    dchSeshOnly.seshName = dch.seshes{i};
    dchSeshOnly.decoherenceTime = dch.decoherenceTime{i};
    dchSeshOnly.decoherenceStartDelay = dch.decoherenceStartDelay{i};
    dchSeshOnly.decoherenceTimeIdx = dch.decoherenceTimeIdx{i};
    dchSeshOnly.Fs = dch.Fs;
    
    seshIndx = ismember(allCells.metadata(:,1),seshes{i});
    nCells = sum(seshIndx);
   
    % Assign struct shaped data to cells
    for m = 1:nCells
       allCells_Dch(z+m).dch = dchSeshOnly;
    end

    z = z + nCells;
    fprintf('Session: %d; Adding dch %d for %d/%d cells\n', i,nCells,z,count)
end

allCells.dch = allCells_Dch;