 function plotAllCells(allCells,paramsPath)
%% Make Plots. Use after running Pool All cells
% Run all the plotting functions
% input: ketamineCells struct

% add plotting functions to path
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath(pwd))

% if ~exist('paramsPath','var')
%     params = readtable('./UniversalParams.xlsx');
% else
%     params = readtable(paramsPath);
% end

%% Filter cells

% Filter for cells in session where only ketamine was delivered
seshIndx = ismember(allCells.metadata(:,8),'ketamine');
ketamineCells = filterAllCellsStruct(allCells,seshIndx);
fprintf('done filtering ketamineCells\n');

% filter for WT mec cells with ketamine
seshIndx = ismember(ketamineCells.metadata(:,4),'WT');
wt_ket_Cells = filterAllCellsStruct(ketamineCells,seshIndx);
fprintf('done filtering for WT-ket cells\n');

% filter for WT mec cells with ketamine not interurons
seshIndx = ~wt_ket_Cells.interneuronFlag;
wt_ket_Cells_noInterneurons = filterAllCellsStruct(wt_ket_Cells,seshIndx);
fprintf('done filtering for WT-ket cells noInterneurons\n');

% filter for WT mec cells with ketamine only interurons
seshIndx = wt_ket_Cells.interneuronFlag;
wt_ket_Cells_onlyInterneurons = filterAllCellsStruct(wt_ket_Cells,seshIndx);
fprintf('done filtering for WT-ket cells only Interneurons\n');

% filter for WT mec cells with gain change cells
seshIndx = logical(wt_ket_Cells_noInterneurons.gainModulationValues(:,4));
wt_ket_Cells_gainChange = filterAllCellsStruct(wt_ket_Cells_noInterneurons,seshIndx);
fprintf('done filtering for WT-ket cells only gain change\n');

% filter for WT mec cells with not gain change cells
seshIndx = logical(~wt_ket_Cells_noInterneurons.gainModulationValues(:,4));
wt_ket_Cells_noGainChange = filterAllCellsStruct(wt_ket_Cells_noInterneurons,seshIndx);
fprintf('done filtering for WT-ket cells only Interneurons\n');


% WT Cells with only stable Cells
wt_ket_onlyStable = filterAllCellsStruct(wt_ket_Cells,wt_ket_Cells.stabilityFlag);

% filter for stable gain change cells that are not interneurons
wt_ket_Cells_noInterneurons_onlyStable = filterAllCellsStruct(wt_ket_Cells_noInterneurons,wt_ket_Cells_noInterneurons.stabilityFlag);
seshIndx = logical(wt_ket_Cells_noInterneurons_onlyStable.gainModulationValues(:,4));
wt_ket_Cells_stableGainChange = filterAllCellsStruct(wt_ket_Cells_noInterneurons_onlyStable,seshIndx);
fprintf('done filtering for WT-ket cells stable only gain change\n');

%% Example session
seshIndx = ismember(allCells.metadata(:,1),'G3_190705_baseline1+controlinjx1+ketamine1_fr+corr');
singleSession_Cells = filterAllCellsStruct(allCells,seshIndx);

%% ========= HCN1KO FILTERING ========
seshIndx = ismember(ketamineCells.metadata(:,4),'KO');
hcn1ko_ket_Cells = filterAllCellsStruct(ketamineCells,seshIndx);
fprintf('done filtering for HCN1ko-ket cells\n');

% filter for HCN1KO mec cells with ketamine not interurons
seshIndx = ~hcn1ko_ket_Cells.interneuronFlag;
hcn1ko_ket_Cells_noInterneurons = filterAllCellsStruct(hcn1ko_ket_Cells,seshIndx);
fprintf('done filtering for hcn1ko-ket cells noInterneurons\n');

% filter for HCN1KO mec cells with ketamine only interurons
seshIndx = hcn1ko_ket_Cells.interneuronFlag;
hcn1ko_ket_Cells_onlyInterneurons = filterAllCellsStruct(hcn1ko_ket_Cells,seshIndx);
fprintf('done filtering for hcn1ko-ket cells only Interneurons\n');

% filter for hcn1ko mec cells with gain change cells
seshIndx = logical(hcn1ko_ket_Cells_noInterneurons.gainModulationValues(:,4));
hcn1ko_ket_Cells_gainChange = filterAllCellsStruct(hcn1ko_ket_Cells_noInterneurons,seshIndx);
fprintf('done filtering for hcn1ko-ket cells only gain change\n');

% filter for hcn1ko mec cells with not gain change cells
seshIndx = logical(~hcn1ko_ket_Cells_noInterneurons.gainModulationValues(:,4));
hcn1ko_ket_Cells_noGainChange = filterAllCellsStruct(hcn1ko_ket_Cells_noInterneurons,seshIndx);
fprintf('done filtering for hcn1ko-ket cells only Interneurons\n');

% filter for hcn1ko stable cells
hcn1ko_ket_Cells_onlyStable = filterAllCellsStruct(hcn1ko_ket_Cells,hcn1ko_ket_Cells.stabilityFlag);
fprintf('done filtering for hcn1ko-ket stable cells\n');

% ========= MK801 & Control FILTERING ========
seshIndx = ismember(allCells.metadata(:,8),'MK801'); 
MK801_cells = filterAllCellsStruct(allCells,seshIndx);
fprintf('done filtering for MK801 cells\n');
% filter for WT mec cells with ketamine
seshIndx = ismember(MK801_cells.metadata(:,4),'WT');
wt_mk801_Cells = filterAllCellsStruct(MK801_cells,seshIndx);
fprintf('done filtering for WT-mk801 cells\n');
% filter for MK801 cells with HCN1ko
seshIndx = ismember(MK801_cells.metadata(:,4),'KO');
hcn1ko_mk801_Cells = filterAllCellsStruct(MK801_cells,seshIndx);
fprintf('done filtering for HCN1ko-mk801 cells\n');

% Control
seshIndx = ismember(allCells.metadata(:,8),'control');
controlCells = filterAllCellsStruct(allCells,seshIndx);
seshIndx = ismember(controlCells.metadata(:,4),'WT');
controlCells = filterAllCellsStruct(controlCells,seshIndx);
fprintf('done filtering Control Cells\n');






% ========  Example session with only Stable Cells======== 
% WT 
%'G3_190704_baseline1+controlinjx1+ketamine1_fr+corr'
% 'G3_190705_baseline1+controlinjx1+ketamine1_fr+corr'
% 'G4_190620_baseline1+controlinjx1+ketamine1_fr+corr'
% 'G5_190708_baseline1+controlinjx1+ketamine1_fr+corr'
% 'HCNd1_190808_baseline1+controlinjx1+ketamine1_fr+corr'

singleSession_Cells_onlyStable = filterAllCellsStruct(singleSession_Cells,singleSession_Cells.stabilityFlag);
wt_ket_Cells_onlyStable = filterAllCellsStruct(wt_ket_Cells,wt_ket_Cells.stabilityFlag);
hcn1ko_ket_Cells_onlyStable = filterAllCellsStruct(hcn1ko_ket_Cells,hcn1ko_ket_Cells.stabilityFlag);

fprintf('=====\ndone filtering\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Single Cell Raster plots with Decoherence Period Highlighted
% by trials and then combine them into combined session rasters
save_figs = true;
image_save_dir = '../fkm_analysis/rasters_dch';
plotSaveCombine_SingleCellRastersPlotsWithDecoherence(allCells,image_save_dir,save_figs)

%% plots All cells rasters & ratemaps
save_figs = true;
image_save_dir = '../fkm_analysis/rasters';
plotAllSingleCells(allCells,image_save_dir,save_figs)

% Combine single cell rasters into large session pngs 
size_vert = 1250; % changed from 1042 to 346
size_horiz = 2083; % changed from 333 to 667 for repeating tracks
combinebySesh(allCells,image_save_dir,size_vert,size_horiz)

%% plots all single cells rasters only
save_figs = true;
image_save_dir = '../fkm_analysis/rasters_black&white';
plotAllSingleCellsRastersOnly(allCells,image_save_dir,save_figs)

% Combine single cell rasters into large session pngs 
size_vert = 1042; 
size_horiz = 333;
combinebySesh(allCells,image_save_dir,size_vert,size_horiz)
%% plots all STABLE single cells rasters only
save_figs = true;
image_save_dir = '../fkm_analysis/rasters_black&white/stableOnly';
plotAllSingleCellsRastersOnly(allCells_onlyStable,image_save_dir,save_figs)

% Combine single cell rasters into large session pngs 
size_vert = 1042; 
size_horiz = 333;
combinebySesh(allCells_onlyStable,image_save_dir,size_vert,size_horiz)
%% plots all single cells RATEMAPS only
save_figs = true;
image_save_dir = '../fkm_analysis/ratemaps';
plotAllSingleCellsRatemapsOnly(allCells,image_save_dir,save_figs)

% Combine single cell rasters into large session pngs 
size_vert = 1042; 
size_horiz = 333;
combinebySesh(allCells,image_save_dir,size_vert,size_horiz)

%% plots all STABLE single cells RATEMAPS only
save_figs = true;
image_save_dir = '../fkm_analysis/ratemaps/stableOnly';
plotAllSingleCellsRatemapsOnly(allCells_onlyStable,image_save_dir,save_figs)

% Combine single cell rasters into large session pngs 
size_vert = 1042; 
size_horiz = 333;
combinebySesh(allCells_onlyStable,image_save_dir,size_vert,size_horiz)

%% plots all umap embeddings
savePath = '../fkm_analysis/umap';
save_figs = true;
plot_UMAPdataEmbedding(allCells,savePath,save_figs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEHAVIOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_BehaviorbySesh(cells,savefigTF)
plot_sessionBehaviorStats(cells)
[lickAccuracyByTrial,timewarpedLickAccuracy,trial,post] = plot_SingleSessionBehavior(cells,savefigTF);

plot_timewarped_lickBehavior(cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot Firing Rate over Time 5 min before injection and 10 min after injection
plot_FRoverTime5minBefore10minafter(wt_ket_Cells)

plot_FRoverTime5minBefore10minafter(wt_ket_Cells_noInterneurons)
plot_FRoverTime5minBefore10minafter(wt_ket_Cells_onlyInterneurons)
plot_FRoverTime5minBefore10minafter(wt_ket_Cells_stableGainChange)
%% Plot Stats comparing Firing Rate over Time 5 min before injection and 5 min after injection
plot_STATS_5minBefore5minafter(wt_ket_Cells)
plot_STATS_5minBefore5minafter(wt_ket_Cells_stableGainChange)
plot_STATS_5minBefore5minafter(wt_ket_Cells_onlyInterneurons)

%% Plot Firing Rate over Time 5min before and 60 min after injection
plot_FRneg5to60minAfterKetamineInjx(wt_ket_Cells,'Ketamine-induced Avg FR Change');
plot_FRneg5to60minAfterKetamineInjx(wt_ket_Cells_stableGainChange,'Ketamine-induced Avg FR - Grid Cells');
plot_FRneg5to60minAfterKetamineInjx(wt_ket_Cells_onlyInterneurons,'Ketamine-induced Avg FR - Interneurons');

%% Plot Firing Rate over Time 5min before and 60 min after injection by mouse
plot_neg5to60minAfterKetamineInjx_byMouse(wt_ket_Cells_stableGainChange)
plot_neg5to60minAfterKetamineInjx_byMouse(wt_ket_Cells_onlyInterneurons)
%% Plot Firing Rate over Trials by Mouse
plot_avgFRbyMouse(wt_ket_Cells,'Average FR by Mouse')
plot_avgFRbyMouse(wt_ket_Cells_stableGainChange,'Average FR by Mouse - Grid Cells')
plot_avgFRbyMouse(wt_ket_Cells_onlyInterneurons,'Average FR by Mouse - Interneurons')
%%
plot_HistfitFRscore(wt_ket_Cells,'WT');
plot_HistfitFRscore(wt_ket_Cells_stableGainChange,'WT');
plot_HistfitFRscore(wt_ket_Cells_onlyInterneurons,'WT');
%% Plot Peakiness Curves over Trials
plot_peakinessCurves(wt_ket_Cells)

%% Plot Time Warped Stability Score Curves
plot_timeWarpedStabilityScoreCurves(wt_ket_Cells)

%% Plot Time Warped Stability Score Curves with only Stable Cells
plot_timeWarpedStabilityScoreCurves(wt_ket_Cells_onlyStable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Nice Single Cell Figure
% i = 1148 %  G4, WT
% i = 300 %  G2, WT
% i = 3213%  HCNd1, WT
% i = 3484 %  HCNd2, KO
% i = 545 % G3, WT
% i = 285 % G2, WT
% i = 297 % G2, WT
% i = 302 % G2, WT
% i = 481 % G3, WT
% i = 258 %G2, WT
% i = 297 % G2,WT
% i = 689
% i = 566; %G3, WT
plot_niceSingleCellFig(allCells,566); % G3_190704_baseline1+controlinjx1+ketamine1_fr+corr
plot_niceSingleCellFig(allCells,568); % G3_190704_baseline1+controlinjx1+ketamine1_fr+corr
plot_niceSingleCellFig(allCells,285); 
%% Plot Correlation Matrix by sessions
% Example session
seshIndx = ismember(allCells.metadata(:,1),'G3_190704_baseline1+controlinjx1+ketamine1_fr+corr');
singleSession_Cells = filterAllCellsStruct(allCells,seshIndx);
plot_correlationMatrixBySessionThenAvg(singleSession_Cells);

% Example session2
seshIndx = ismember(allCells.metadata(:,1),'G4_190622_baseline1+controlinjx1+ketamine1_fr+corr');
singleSession_Cells = filterAllCellsStruct(allCells,seshIndx);
plot_correlationMatrixBySessionThenAvg(singleSession_Cells);

% Example session3
seshIndx = ismember(allCells.metadata(:,1),'G1_190817_baseline1+controlinjx1+ketamine1_fr+corr');
singleSession_Cells = filterAllCellsStruct(allCells,seshIndx);
plot_correlationMatrixBySessionThenAvg(singleSession_Cells);

% PCA for example session G3_190704_baseline1+controlinjx1+ketamine1_fr+corr
plot_PCA(singleSession_Cells, false);

% Plot Correlation Matrix
plot_correlationMatrix(singleSession_Cells_onlyStable, 'WT');
plot_correlationMatrix(singleSession_Cells, 'WT');

% 
plot_correlationScoreCurves(wt_ket_Cells_onlyStable,'WT')
plot_correlationScoreCurves(wt_ket_Cells,'WT')



%% PLOT Histfit on Corr score
plot_HistfitKetCorrEffectScore(wt_ket_Cells, 'WT')
plot_HistfitKetCorrEffectScore(wt_ket_Cells_onlyStable, 'WT')
plot_HistfitKetCorrEffectScore(wt_ket_Cells_noInterneurons,'WT')
plot_HistfitKetCorrEffectScore(wt_ket_Cells_onlyInterneurons,'WT')
%% Plot Correlation Matrix by sessions
plot_correlationMatrixBySessionThenAvg(wt_ket_Cells);


%% Decoherence Plots
plotDecoherencePlots(wt_ket_Cells);

%% Plot PCA by session for WT animals
plot_PCAbySesh(allCells,false);

%% Plot FR comparisons of decoherence 
plot_cells1_vs_cells2_dchPeriod_STATS(wt_ket_Cells_gainChange,'Grid cells',wt_ket_Cells_onlyInterneurons,'interneurons')
plot_cells1_vs_cells2_dchPeriod_STATS(wt_ket_Cells_gainChange,'Gain Change',wt_ket_Cells_noGainChange,'No Gain Change')

%% Comparing Interneurons FR time course 
plot_FRoverTime5minBefore10minafter(wt_ket_Cells_onlyInterneurons)
plot_FRoverTime5minBefore10minafter(wt_ket_Cells_noInterneurons)

plot_FRoverTime5minBefore10minafter(wt_ket_Cells_gainChange)
plot_FRoverTime5minBefore10minafter(wt_ket_Cells_noGainChange)
%%
plot_STATS_cells1_vs_cells2_frDiff_5minBefore5minafter(wt_ket_Cells_gainChange,'Grid cells',wt_ket_Cells_onlyInterneurons,'interneurons')
plot_STATS_cells1_vs_cells2_frDiff_5minBefore5minafter(wt_ket_Cells_gainChange,'Gain Change',wt_ket_Cells_noGainChange,'No Gain Change')

%%
plot_STATS_5minBefore5minafter(wt_ket_Cells_onlyInterneurons)
plot_STATS_5minBefore5minafter(wt_ket_Cells_noInterneurons)

plot_STATS_5minBefore5minafter(wt_ket_Cells_gainChange)
plot_STATS_5minBefore5minafter(wt_ket_Cells_noGainChange)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Attractor Network Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

findAttractorDynamicsbySession(wt_ket_Cells_gainChange,true,"wt_cell_gainchange_only")
findAttractorDynamicsbySession(wt_ket_Cells_noInterneurons,true,"ket_noInterneurons")
findAttractorDynamicsbySession(wt_ket_Cells_stableGainChange,true,"wt_cell_gainchange_only_stableOnly")
findAttractorDynamicsbySession(wt_ket_onlyStable,true,"wt_cell_stableOnly")
findAttractorDynamicsbySession(wt_ket_Cells_onlyInterneurons,true,"wt_cell_ket_interneuronOnly")

findAttractorDynamicsbySession(hcn1ko,true,"hcn1ko_ket_Cells_onlyStable")
findAttractorDynamicsbySession(wt_mk801_Cells,true,"wt_mk801_Cells")
findAttractorDynamicsbySession(hcn1ko_mk801_Cells,true,"hcn1ko_mk801_Cells")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE HCN1ko
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Firing Rate over Time 5 min before injection and 10 min after injection
plot_FRoverTime5minBefore10minafter(hcn1ko_ket_Cells)

%% Plot Stats comparing Firing Rate over Time 5 min before injection and 5 min after injection
plot_STATS_5minBefore5minafter(hcn1ko_ket_Cells)

%% Plot Peakiness Curves over Trials
plot_peakinessCurves(hcn1ko_ket_Cells_onlyStable)

%% Plot Time Warped Stability Score Curves
plot_timeWarpedStabilityScoreCurves(hcn1ko_ket_Cells)

%% Plot Time Warped Stability Score Curves with only Stable Cells
plot_timeWarpedStabilityScoreCurves(hcn1ko_ket_Cells_onlyStable)

%% Plot Correlation Score Curves
plot_correlationScoreCurves(wt_ket_Cells_onlyStable,'WT')
plot_correlationScoreCurves(hcn1ko_ket_Cells_onlyStable,'HCN1ko')
plot_correlationScoreCurveComparison(wt_ket_Cells_onlyStable,hcn1ko_ket_Cells_onlyStable);
% plot_correlationScoreCurves(ketamineCells, KOcellsIndx,'KO')
%% Plot Timewarped Correlation Score Curves
plot_timeWarpedCorrelationScoreCurves(wt_ket_Cells_onlyStable);
plot_timeWarpedCorrelationScoreCurves(hcn1ko_ket_Cells_onlyStable);

%% Plot Nice Single Cell Figure
plot_niceSingleCellFig(hcn1ko_ket_Cells,[]);
%% Plot Firing Rate over Time 5 min before injection and 10 min after injection
plot_FRoverTime5minBefore10minafter(hcn1ko_ket_Cells);
%% Plot Stats comparing Firing Rate over Time 5 min before injection and 5 min after injection
% plot_STATS_5minBefore5minafter(hcn1ko_ket_Cells)

%% Plot Firing Rate over Time 5min before and 60 min after injection
plot_FRneg5to60minAfterKetamineInjx(hcn1ko_ket_Cells,'Ketamine-induced Avg FR Change on HCN1ko');
 plot_FRneg5to60minAfterKetamineInjx(hcn1ko_ket_Cells_onlyInterneurons,'Ketamine-induced Avg FR Change - interneurons HCN1ko');
  plot_FRneg5to60minAfterKetamineInjx(hcn1ko_ket_Cells_gainChange,'Ketamine-induced Avg FR Change - Grid Cells HCN1ko');
%% Plot Firing Rate over Trials by Mouse
plot_avgFRbyMouse(hcn1ko_ket_Cells);
%% Plot Firing Rate over Time 5min before and 60 min after injection by mouse
plot_neg5to60minAfterKetamineInjx_byMouse(hcn1ko_ket_Cells)
%% Plot dch plot
plotDecoherencePlots(hcn1ko_ket_Cells);

%% 
plot_dch_autocorrelationScore_cells1_VS_cells2(wt_ket_Cells,hcn1ko_ket_Cells, 'WT Ketamine','HCN1ko Ketamine')
plot_dch_autocorrelationScore_cells1_VS_cells2(wt_ket_Cells,wt_mk801_Cells, 'WT Ketamine','MK801 Ketamine')
plot_dch_autocorrelationScore_cells1_VS_cells2(hcn1ko_ket_Cells,wt_mk801_Cells, 'HCN1ko Ketamine','MK801 Ketamine')
plot_dch_autocorrelationScore_cells1_VS_cells2(hcn1ko_mk801_Cells,wt_mk801_Cells, 'HCN1ko-MK801 Ketamine','MK801 Ketamine')

%% Plot Correlation Score Curve Comparisions
plot_correlationScoreCurveComparison(wt_ket_onlyStable,hcn1ko_ket_Cells_onlyStable)

%% Plot Timewarped Correlation Score Curve Comparisions
plot_correlationScoreCurveComparison(wt_ket_onlyStable,hcn1ko_ket_Cells_onlyStable)

%% Plot Peakiness Curves over Trials
plot_peakinessCurves(hcn1ko_ket_Cells_onlyStable)

%% PLOT Histfit on Corr score
plot_HistfitKetCorrEffectScore(hcn1ko_ket_Cells, 'KO')


%% Plot Correlation Matrix
plot_correlationMatrix(wt_ket_Cells, 'WT');
plot_correlationMatrix(wt_ket_Cells_onlyInterneurons, 'WT');
plot_correlationMatrix(wt_ket_Cells_noInterneurons, 'WT');
plot_correlationMatrix(wt_ket_Cells_gainChange, 'WT');
plot_correlationMatrix(wt_ket_Cells_noGainChange, 'WT');

plot_correlationMatrix(hcn1ko_ket_Cells, 'KO');
plot_correlationMatrix(hcn1ko_ket_Cells_onlyInterneurons, 'KO');
plot_correlationMatrix(hcn1ko_ket_Cells_noInterneurons, 'KO');
plot_correlationMatrix(hcn1ko_ket_Cells_gainChange, 'KO');
plot_correlationMatrix(hcn1ko_ket_Cells_noGainChange, 'KO');


% 
% hcn1ko_ket_Cells_noGainChange

plot_correlationMatrixBySessionThenAvg(hcn1ko_ket_Cells);
plot_correlationMatrixBySessionThenAvg(hcn1ko_ket_Cells_gainChange);
plot_correlationMatrixBySessionThenAvg(hcn1ko_ket_Cells_noGainChange);
plot_correlationMatrixBySessionThenAvg(hcn1ko_ket_Cells_onlyInterneurons);
plot_correlationMatrixBySessionThenAvg(hcn1ko_ket_Cells_noInterneurons);
%% Plot Behavior by Sessions for KO ketamine animals
plot_BehaviorbySesh(fltrCells,true);


%% Plot FR comparisons of decoherence 
plot_cells1_vs_cells2_dchPeriod_STATS(hcn1ko_ket_Cells_gainChange,'HCN1ko Grid cells',hcn1ko_ket_Cells_onlyInterneurons,'HCN1ko interneurons')
plot_cells1_vs_cells2_dchPeriod_STATS(hcn1ko_ket_Cells_gainChange,'Gain Change',hcn1ko_ket_Cells_noGainChange,'No Gain Change')

%% Comparing Interneurons FR time course 
plot_FRoverTime5minBefore10minafter(hcn1ko_ket_Cells_onlyInterneurons)
plot_FRoverTime5minBefore10minafter(hcn1ko_ket_Cells_noInterneurons)

plot_FRoverTime5minBefore10minafter(hcn1ko_ket_Cells_gainChange)
plot_FRoverTime5minBefore10minafter(hcn1ko_ket_Cells_noGainChange)
%%
plot_STATS_cells1_vs_cells2_frDiff_5minBefore5minafter(hcn1ko_ket_Cells_noInterneurons,'HCN1ko cells',hcn1ko_ket_Cells_onlyInterneurons,'interneurons')
plot_STATS_cells1_vs_cells2_frDiff_5minBefore5minafter(hcn1ko_ket_Cells_gainChange,'Gain Change',hcn1ko_ket_Cells_onlyInterneurons,'interneurons')
plot_STATS_cells1_vs_cells2_frDiff_5minBefore5minafter(hcn1ko_ket_Cells_gainChange,'Gain Change',hcn1ko_ket_Cells_noGainChange,'No Gain Change')
%%
plot_STATS_5minBefore5minafter(hcn1ko_ket_Cells_onlyInterneurons)
plot_STATS_5minBefore5minafter(hcn1ko_ket_Cells_noInterneurons)

plot_STATS_5minBefore5minafter(hcn1ko_ket_Cells_gainChange)
plot_STATS_5minBefore5minafter(hcn1ko_ket_Cells_noGainChange)

%%
plot_STATS_cells1_vs_cells2_frDiff_5minBefore5minafter(wt_ket_Cells_gainChange,'WT Grid cells',hcn1ko_ket_Cells_gainChange,'HCN1ko Grid Cells')
plot_cells1_vs_cells2_dchPeriod_STATS(wt_ket_Cells_gainChange,'WT Grid cells',hcn1ko_ket_Cells_gainChange,'HCN1ko Grid Cells')
plot_cells1_vs_cells2_dchPeriod_STATS(wt_ket_Cells_onlyInterneurons,'WT Interneurons',hcn1ko_ket_Cells_onlyInterneurons,'HCN1ko Interneurons')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot MK801 WT Cells

%% Plot Raster Grid
plotRasterGrid(fltrCells,100)
% Plot Fr Grid
plotFRGrid(fltrCells,100)

%% Plot Nice Single Cell Figure
plot_niceSingleCellFig(fltrCells);

%% Plot Firing Rate over Time 5min before and 60 min after injection
plot_FRneg5to60minAfterKetamineInjx(fltrCells);

%% Plot Firing Rate over Trials by Mouse
plot_avgFRbyMouse(fltrCells)

%% Plot Correlation Score Curves
plot_correlationScoreCurves(fltrCells,'WT')
% plot_correlationScoreCurves(ketamineCells, KOcellsIndx,'KO')

%% Plot Correlation Score Curve Comparisions
plot_correlationScoreCurveComparison(fltrCells, ko_fltrCells)


%% Plot Peakiness Curves over Trials
plot_peakinessCurves(fltrCells)

%% PLOT Histfit on Corr score
plot_HistfitKetCorrEffectScore(fltrCells, 'WT-MK801')


%% Plot Correlation Matrix
plot_correlationMatrix(fltrCells, 'MK801');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MISC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Behavior by Sessions for WT animals
plot_BehaviorbySesh(wt_ket_Cells,true);

%% Plot Behavior WT vs Ketamine


%% Filter cells to ketamine cells
seshIndx = ismember(allCells.metadata(:,8),'control');
fltrCells = filterAllCellsStruct(allCells,seshIndx);
fprintf('done filtering ketamineCells\n');


%% Plot peakiness sessions
seshes = unique(cellfun(@num2str,fltrCells.metadata(:,1),'uni',0));

rows = ceil(sqrt(numel(seshes)));
cols = rows;
figure(); hold on;
for i = 1:numel(seshes)
    
    seshIndx = ismember(fltrCells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(fltrCells,seshIndx);
    handles.(sprintf('sp%d',i)) = subplot(rows,cols,i);
    try
        plot_peakinessSubset(seshCells, handles.(sprintf('sp%d',i)))
    catch
        warning('Could not do %d',i);
        continue
    end
end

%% Plot Stability Table STats
plot_stabilityTableStats(ketamineCells, stabilityTable);

%% PLOT Histfit on FR score
plot_HistfitFRscore(ketamineCells, WTcellsIndx,'WT');
plot_HistfitFRscore(ketamineCells, KOcellsIndx,'KO');
 
plot_HistfitFRscore(ketamineCells, WTtotalStabilityIndx,'WT');
plot_HistfitFRscore(ketamineCells, KOtotalStabilityIndx,'KO');
%% Plot Distribution of Ketamine Correlation Scores
plot_HistfitKetCorrEffectScore(ketamineCells, WTcellsIndx,'WT')
plot_HistfitKetCorrEffectScore(ketamineCells, KOcellsIndx,'KO')

%% Plot Stability Score Curve
plot_stabilityScore(ketamineCells, [], 'all cells')
plot_stabilityScore(ketamineCells, WTcellsIndx,'WT')
plot_stabilityScore(ketamineCells, KOcellsIndx,'KO')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRELATION MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = fieldnames(allSpatialIndx);

sessionMap = [];
count = 0;
for k=1:numel(fn)
    if(isnumeric(allSpatialIndx.(fn{k})))
        seshCellNum = size(allSpatialIndx.(fn{k}),2);
        sessionMap(count+1:count+seshCellNum,1) = k;
        count = count + seshCellNum;
    end
end
%%
for k=1:numel(fn)
    try
    testCells = ketamineCells.spatialFR(sessionMap == k,:,:);

    
    numCells = size(testCells,1);
    trialNum = size(testCells,2);
    spatialBins = size(testCells,3);

    flatFR2 = reshape(permute(testCells,[3 1 2]), [numCells*spatialBins, trialNum]);

    %P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
    corrMatrix = corr(fillmissing(flatFR2,'linear')); 
    figure(1); clf;
    imagesc(corrMatrix); colorbar;
    title(fn{k},'Interpreter', 'none');
    set(gca,'TickDir','out');
    set(gca,'ticklength',[0.005 0.025]);
    set(gca,'layer','bottom');
    box off;
    axis square;
    set(gca,'FontSize',30);
    set(gca,'FontName','Helvetica');
%     set(gcf,'Position',[100 100 1000 1000])

    tempAllCellScorrMatrix = nan(numCells, trialNum,trialNum);
    for i = 1:size(testCells,1)
        %calculate trial by trial correlation matrix for one cell 
        singleCellFR = squeeze(testCells(i,:,:))';
        corrMatrix = corr(fillmissing(singleCellFR,'linear'));
        tempAllCellScorrMatrix(i, :, :) = corrMatrix;
%         imagesc(corrMatrix)
%         pause 
    end
    figure(2); clf;
    imagesc(squeeze(nanmean(tempAllCellScorrMatrix,1))); colorbar;
    title(fn{k},'Interpreter', 'none');
    set(gca,'TickDir','out');
    set(gca,'ticklength',[0.005 0.025]);
    set(gca,'layer','bottom');
    box off;
    axis square;
    set(gca,'FontSize',30);
    set(gca,'FontName','Helvetica');
%     set(gcf,'Position',[100 100 1000 1000])
   pause 
    end
end

%%

testCells = ketamineCells.spatialFR10;


numCells = size(testCells,1);
trialNum = size(testCells,2);
spatialBins = size(testCells,3);

flatFR2 = reshape(permute(testCells,[3 1 2]), [numCells*spatialBins, trialNum]);

%P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
% corrMatrix = corr(fillmissing(flatFR2,'linear')); 
corrMatrix = corr(ketamineCells.spatialFR10); 

figure(1); clf;
imagesc(corrMatrix); colorbar;
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box on;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Trial by Trial Population Activity Correlation Matrix(%s)',filter))


figure(2); clf;
imagesc(squeeze(nanmean(ketamineCellsCorrMatrix,1)),[0, 0.2]); colorbar; 
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box on;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Avg Cell Trial by Trial Correlation Matrix(%s)',filter))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MISC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Sort by peakiness
% [~,sortIndx] = sortrows(...
%     nanmean(ketamineCells.peakiness,2)/max(nanmean(ketamineCells.peakiness,2))...
%     + (nanmean(nanmean(ketamineCells.spatialFR10,2),3))/max(nanmean(nanmean(ketamineCells.spatialFR10,2),3))...
%     ,'descend');
[~,sortIndx] = sortrows(nanmean(wt_ket_Cells.peakiness,2),'ascend');
sortedCells = sortAllCellsStruct(wt_ket_Cells,sortIndx);
%% Plots single cells sorted by peakiness
image_save_dir = '../fkm_analysis/rasters_by_peakiness';
plotAllSingleCells(sortedCells,image_save_dir,false)

%% Plot Raster Grid
plotRasterGrid(wt_ket_Cells,100)
% Plot Fr Grid
plotFRGrid(wt_ket_Cells,100)

%%
plot_sortedMatrix(wt_ket_Cells.peakiness,sortIndx,'ascend')
%% Sort by Peakiness and Firing Rate - Aggregated
sortIndx = sortByTwoCols(wt_ket_Cells);
sortedCells = sortAllCellsStruct(wt_ket_Cells,sortIndx);
plotRasterGrid(wt_ket_Cells,36)

%% Generate Indices
% 
% WTcellsIndx = strcmp(ketamineCells.metadata(:,4), 'WT')';
% KOcellsIndx = strcmp(ketamineCells.metadata(:,4), 'KO')';

% stabilityTable = findStableCells(allCells); % {'totalStability', 'baselineStability', 'acuteDrugStability', 'endingStability', 'gainStability'}
% stabilityThreshold = 0.2;
% totalStabilityIndx = stabilityTable.totalStability > stabilityThreshold;

% stabilityTable = findStableCells(ketamineCells);
% stabilityThreshold = 0.2;
% WTtotalStabilityIndx = (stabilityTable.totalStability > stabilityThreshold) & WTcellsIndx;
% KOtotalStabilityIndx = (stabilityTable.totalStability > stabilityThreshold) & KOcellsIndx;

% baselineStabilityIndx = stabilityTable.baselineStability > stabilityThreshold;
% WTbaselineStabilityIndx = (stabilityTable.baselineStability > stabilityThreshold) & WTcellsIndx;
% KObaselineStabilityIndx = (stabilityTable.baselineStability > stabilityThreshold) & KOcellsIndx;
end
