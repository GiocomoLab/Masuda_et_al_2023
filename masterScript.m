% Masters Script for Ketamine Paper
% Francis Kei Masuda 2020

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
sessionMetaDataPath = '/Users/KeiMasuda/Desktop/fkm_analysis/SessionList.xlsx';
allCells = poolAllCells(filter,sessionMetaDataPath);

% % Filter for cells in session where only ketamine was delivered
% seshIndx = ismember(allCells.metadata(:,8),'ketamine');
% ketamineCells = filterAllCellsStruct(allCells,seshIndx);
% fprintf('done filtering ketamineCells\n');
% 
% % Filter for only WT mec cells with ketamine
% seshIndx = ismember(ketamineCells.metadata(:,4),'WT');
% fltrCells = filterAllCellsStruct(ketamineCells,seshIndx);
% fprintf('done filtering for WT cells\n');

%% Calculate a decoherence band for each session in given cells
dch = calcDecoherenceBandForSessions(allCells);
save('/Users/keimasuda/Desktop/fkm_analysis/dch.mat','dch')
%%
load('/Users/keimasuda/Desktop/fkm_analysis/dch.mat')

%% Plot Single Cell Raster plots with Decoherence Period Highlighted with
% Trial-based Dch Index and then combine them into combined session rasters
save_figs = true;
image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/rasters_dch';
plotSaveCombine_SingleCellRastersPlotsWithDecoherence(cells,dch,image_save_dir,save_figs)

%%
plotAllCells(allCells,dch);