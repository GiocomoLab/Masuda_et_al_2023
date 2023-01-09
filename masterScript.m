% Master Script for Ketamine Paper
% Francis Kei Masuda 2020
% Requires JohnKeiNPAnalysis to be on the same directory level as fkm_analysis

addpath('./plottingFxns')
addpath(genpath(pwd))
paramsPath = './UniversalParams.xlsx';
filter = 'mec';
%% Only run once after each session has been spike sorted
% run stitch sessions together, plot raster plots, combine rasters, plot
% lick data for 3 unity sessions. NEED TO UPDATE WITH NEW SYNC DRIFT
% CORRECT

% runMultiPostProcessing
%%
% combinedSessionsPath = '../fkm_analysis/combinedSesh/*.mat';
% saveDir = '../fkm_analysis/combinedSesh/fr_data_matrices_noSmoothing';
% 
% runMultiAnalysis(filter,combinedSessionsPath,saveDir,paramsPath);

%% Generate spatial indx to figure out what cells exist to pool together in the next step
% Not needed if a spatial index has already be created
spatialIndx = generateSpatialIndx(filter);
%% Pool All the Cells from calculated metadata files from sessions identified in spreadsheet into one big struct
tic
sessionMetaDataPath = '../fkm_analysis/SessionList.xlsx';
allCells = poolAllCells(filter,sessionMetaDataPath);
toc
%% Calculate a decoherence band for each session in given cells
% dchFolderPath = '../fkm_analysis/dch/';
% dch = calcDecoherenceBandForSessions(allCells,dchFolderPath,paramsPath);
% dchFilePath = '../fkm_analysis/dch.mat';
% save(dchFilePath,'dch');
%% if dechorence bands have already been calculated — use this to add dch 
% band to the allCells struct
dchFilePath = '../fkm_analysis/dch.mat';
load(dchFilePath);
allCells = add_Dch_to_allCells(allCells,dch);

% Add a stability flag to all the cells
stabilityTable = findStableCells(allCells); % {'totalStability', 'baselineStability', 'acuteDrugStability', 'endingStability', 'gainStability'}
stabilityThreshold = 0.2;
baselineStabilityIndx = stabilityTable.baselineStability > stabilityThreshold;
allCells.stabilityFlag = baselineStabilityIndx;
fprintf('Done adding a stability flag to AllCells\n');

% Add gain change modulation values (comparing last 10 trials)
gainModulationValues = findGainModulatedCells(allCells);
allCells.gainModulationValues = gainModulationValues;
fprintf('Done adding gainModulationValues to AllCells\n');

% Add interneuron flag (baseline mean FR > 15hz)
allCells.interneuronFlag = findInterneurons(allCells);
fprintf('Done adding interneuron flag to AllCells\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-Figures
% Plot decoherence bands stats
plotDecoherenceBandStats(dch)
% 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code to generate paper figures
plotAllCells(allCells,paramsPath);