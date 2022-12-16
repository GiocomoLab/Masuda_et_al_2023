function allCells = add_Dch_to_allCells(allCells,dch)
% Find unique sessions in cells struct
clear dchSeshOnly
seshes = unique(cellfun(@num2str,allCells.metadata(:,1),'uni',0));
count = numel(allCells.metadata(:,1));

clear allCells_Dch
allCells_Dch(count) = struct(); 
z = 0;


for i = 1:numel(seshes)
    
    dchSeshOnly.idxCellArray = dch.idxCellArray{i};
    dchSeshOnly.idxClusterArray = dch.idxClusterArray{i};
    dchSeshOnly.decoherenceIdx = dch.decoherenceIdx{i};
    dchSeshOnly.seshName = dch.seshes{i};
    dchSeshOnly.decoherenceTime = dch.decoherenceTime{i};
    dchSeshOnly.decoherenceStartDelay = dch.decoherenceStartDelay{i};
    dchSeshOnly.decoherenceTimeIdx = dch.decoherenceTimeIdx{i};
    dchSeshOnly.Fs = dch.Fs;
    dchSeshOnly.umapOutput = dch.umapOutput{i};
    try
        dchSeshOnly.timeDownSample = dch.timeDownSample{i};
    catch
        dchSeshOnly.timeDownSample = 100;
    end
    
    seshIndx = ismember(allCells.metadata(:,1),seshes{i});
    nCells = sum(seshIndx);
   
    % Assign struct shaped data to cells
    for m = 1:nCells
       allCells_Dch(z+m).dch = dchSeshOnly;
    end

    z = z + nCells;
%     fprintf('Session: %d; Adding dch %d for %d/%d cells\n', i,nCells,z,count)
end
fprintf('Done with the Dch to all Cells addition\n');
allCells.dch = allCells_Dch;