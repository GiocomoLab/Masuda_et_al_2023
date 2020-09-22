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
combinedSessionsPath = '/Users/KeiMasuda/Desktop/fkm_analysis/combinedSesh/*.mat';
saveDir = '/Users/KeiMasuda/Desktop/fkm_analysis/combinedSesh/fr_data_matrices_noSmoothing';
paramsPath = '/Users/keimasuda/Desktop/JohnKeiNPAnalysis/UniversalParams.xlsx';
runMultiAnalysis(filter,combinedSessionsPath,saveDir,paramsPath);

%% Generate spatial indx to figure out what cells exist to pool together in the next step
% Not needed if a spatial index has already be created
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
dchFilePath = '/Users/keimasuda/Desktop/fkm_analysis/dch.mat';
save(dchFilePath,'dch');
%% if dechorence bands have already been calculated — use this to add dch 
% band to the allCells struct
dchFilePath = '/Users/keimasuda/Desktop/fkm_analysis/dch.mat';
load(dchFilePath);
allCells = add_Dch_to_allCells(allCells,dch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-Figures
% Plot decoherence bands stats
plotDecoherenceBandStats(dch)
%% Plot Single Cell Raster plots with Decoherence Period Highlighted
% by trials and then combine them into combined session rasters
save_figs = true;
image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/rasters_dch';
plotSaveCombine_SingleCellRastersPlotsWithDecoherence(allCells,image_save_dir,save_figs)


%% plots All cells rasters & ratemaps
save_figs = true;
image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/rasters';
plotAllSingleCells(allCells,image_save_dir,save_figs)

% Combine single cell rasters into large session pngs 
size_vert = 1250; % changed from 1042 to 346
size_horiz = 2083; % changed from 333 to 667 for repeating tracks
combinebySesh(allCells,image_save_dir,size_vert,size_horiz)

%% plots all single cells rasters only
save_figs = true;
image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/rasters_black&white';
plotAllSingleCellsRastersOnly(allCells,image_save_dir,save_figs)

% Combine single cell rasters into large session pngs 
size_vert = 1042; 
size_horiz = 333;
combinebySesh(allCells,image_save_dir,size_vert,size_horiz)

%% plots all single cells RATEMAPS only
save_figs = true;
image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/ratemaps';
plotAllSingleCellsRastersOnly(allCells,image_save_dir,save_figs)

% Combine single cell rasters into large session pngs 
size_vert = 1042; 
size_horiz = 333;
combinebySesh(allCells,image_save_dir,size_vert,size_horiz)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code to generate paper figures
plotAllCells(allCells);