%%
% Only run once after each session has been spike sorted
runMultiPostProcessing

%%

% Reads sessions 

filter = 'mec';
runMultiAnalysis

spatialIndx = generateSpatialIndx(filter);
allCells = poolAllCells(filter);