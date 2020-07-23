function plot_multiDimensionalCellArray(A)
    
    nCol   = cellfun('size', A, 1);
    nRow   = cellfun('size', A, 2);
    Output1 = NaN(max(nCol), numel(A));
    Output2 = NaN(numel(A), max(nRow));

    for iA = 1:numel(A)
        try
          idx = A{iA,:};
          Output1(1:size(idx,1), iA) = idx;
          result = Output1;
        catch
          idx = A{iA,:};
          Output2(iA,1:size(idx,2)) = idx;
          result = Output2';
        end
    end
    figure;
    imagesc(result)
end