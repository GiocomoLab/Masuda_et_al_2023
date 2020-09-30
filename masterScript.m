% Master Script for Ketamine Paper
% Francis Kei Masuda 2020

addpath('./plottingFxns')
filter = 'mec';
%% Only run once after each session has been spike sorted
% run stitch sessions together, plot raster plots, combine rasters, plot
% lick data for 3 unity sessions. NEED TO UPDATE WITH NEW SYNC DRIFT
% CORRECT

runMultiPostProcessing
%%
combinedSessionsPath = '../fkm_analysis/combinedSesh/*.mat';
saveDir = '../fkm_analysis/combinedSesh/fr_data_matrices_noSmoothing';
paramsPath = './UniversalParams.xlsx';
runMultiAnalysis(filter,combinedSessionsPath,saveDir,paramsPath);

%% Generate spatial indx to figure out what cells exist to pool together in the next step
% Not needed if a spatial index has already be created
spatialIndx = generateSpatialIndx(filter);
%% Pool All the Cells from calculated metadata files from sessions identified in spreadsheet into one big struct
sessionMetaDataPath = '../fkm_analysis/SessionList.xlsx';
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
dchFolderPath = '../fkm_analysis/dch/';
dch = calcDecoherenceBandForSessions(allCells,dchFolderPath);
dchFilePath = '../fkm_analysis/dch.mat';
save(dchFilePath,'dch');
%% if dechorence bands have already been calculated — use this to add dch 
% band to the allCells struct
dchFilePath = '../fkm_analysis/dch.mat';
load(dchFilePath);
allCells = add_Dch_to_allCells(allCells,dch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-Figures
% Plot decoherence bands stats
plotDecoherenceBandStats(dch)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code to generate paper figures
plotAllCells(allCells,paramsPath);