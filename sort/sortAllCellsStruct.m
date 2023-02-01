function  sortedCells = sortAllCellsStruct(allCells,sortIndx)
    sortedCells = allCells;
    fields = fieldnames(allCells);
    indxSize = numel(sortIndx);
    for i = 1:numel(fieldnames(allCells))
        fprintf(sprintf('Sorted Field: %d/%d\n',i,numel(fieldnames(allCells))))
        allCellsFields = allCells.(fields{i});
        for j = 1:ndims(allCellsFields)
            if size(allCellsFields,j) == indxSize
                if j == 1
                    if ndims(allCellsFields) == 2
                        sortedCells.(fields{i}) = allCellsFields(sortIndx,:);
                    elseif ndims(allCellsFields) == 3
                        sortedCells.(fields{i}) = allCellsFields(sortIndx,:,:);
                    else
                        error('allCells Struct has a field that is the wrong shape to be filtered1:%s',fields{i});
                    end
                    break;
                elseif j == 2
                    sortedCells.(fields{i}) = allCellsFields(:,sortIndx);
                else
                    error('allCells Struct has a field that is the wrong shape to be filtered1:%s',fields{i});
                end
            end
        end
    end
end

