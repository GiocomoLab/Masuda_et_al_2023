function data = extractSessionValueFromCellsStruct(dataStruct)

data = squeeze(cell2mat(struct2cell(dataStruct)));
data = data(:,1);