function data = extractSessionValueFromCellsStruct(dataStruct)
% gets the behavioral values by extracting the first cell's behavior data
% assumes that the input is all cells from the same session
data = squeeze(cell2mat(struct2cell(dataStruct)));
data = data(:,1);