function result = convertMultidimensionalCellArray2paddedMatrix(A)
% takes a n x 1 cell array where n can be different sizes of m x 1 cell
% arrays and converts it to a matrix padded with nans 

nCol   = cellfun('size', A, 1);

result = NaN(max(nCol), numel(A));

for iA = 1:numel(A)
  idx = A{iA,:};
  result(1:size(idx,1), iA) = idx;
end

