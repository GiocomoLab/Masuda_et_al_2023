function  fltrFields = filterAllCellsStruct(allCells,indx)
    fields = fieldnames(allCells);
    indxSize = numel(indx);
    for i = 1:numel(fieldnames(allCells))
%         fprintf(sprintf('\nFiltering: %s',fields{i}));
        allCellsFields = allCells.(fields{i});
        for j = 1:ndims(allCellsFields)
            if size(allCellsFields,j) == indxSize
                if j == 1
                    if ndims(allCellsFields) == 2
                        fltrFields.(fields{i}) = allCellsFields(indx,:);
                    elseif ndims(allCellsFields) == 3
                        fltrFields.(fields{i}) = allCellsFields(indx,:,:);
                    else
                        error('allCells Struct has a field that is the wrong shape to be filtered1:%s',fields{i});
                    end
                    break;
                elseif j == 2
                    fltrFields.(fields{i}) = allCellsFields(:,indx);
                else
                    error('allCells Struct has a field that is the wrong shape to be filtered1:%s',fields{i});
                end
            end
        end
    end
end