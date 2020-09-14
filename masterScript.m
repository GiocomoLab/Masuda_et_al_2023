%%

addpath('./plottingFxns')
%% Only run once after each session has been spike sorted
% run stitch sessions together, plot raster plots, combine rasters, plot
% lick data for 3 unity sessions. NEED TO UPDATE WITH NEW SYNC DRIFT
% CORRECT

runMultiPostProcessing
%%
filter = 'mec';
%%
combinedSessionsPath = '/Users/KeiMasuda/Desktop/fkm_analysis/combinedSesh/*.mat';
runMultiAnalysis(filter,combinedSessionsPath);

%%
spatialIndx = generateSpatialIndx(filter);
%% Pool All the Cells from calculated metadata files from sessions identified in spreadsheet into one big struct
tic
sessionMetaDataPath = '/Users/KeiMasuda/Desktop/fkm_analysis/SessionList.xlsx';
allCells = poolAllCells(filter,sessionMetaDataPath);

% Filter for cells in session where only ketamine was delivered
seshIndx = ismember(allCells.metadata(:,8),'ketamine');
ketamineCells = filterAllCellsStruct(allCells,seshIndx);
fprintf('done filtering ketamineCells\n');

% Filter for only WT mec cells with ketamine
seshIndx = ismember(ketamineCells.metadata(:,4),'WT');
fltrCells = filterAllCellsStruct(ketamineCells,seshIndx);
fprintf('done filtering for WT cells\n');
toc
%% Calculate a decoherence band for each session in given cells

dch = calcDecoherenceBandForSessions(fltrCells);

%%
plotDecoherenceBandStats(dch)

 %% Plot Single Cell Raster plots with Decoherence Period Highlighted with Trial-based Dch Index
save_figs = true;
image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/rasters_dch';
for i = 1:numel(seshes)
    close all;
    seshIndx = ismember(fltrCells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(fltrCells,seshIndx);
    
    plotDchRaster(decoherenceIdx,seshCells,i, save_figs,image_save_dir)
end
%% Combine single cell rasters into large session pngs 
combinebySesh(fltrCells,image_save_dir)
%%  Plot Single Cell Raster plots with Decoherence Period Highlighted with Spatially-binned Dch Index
% save_figs = false;
% for i = 1:numel(seshes)
%     seshIndx = ismember(fltrCells.metadata(:,1),seshes{i});
%     seshCells = filterAllCellsStruct(fltrCells,seshIndx);
%     plotDchRaster_timeBinnedFR(decoherenceIdx,cells,seshID, Fs, save_figs)
% end

%%
plotAllCells(allCells);