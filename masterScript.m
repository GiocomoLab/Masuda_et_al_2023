%%

addpath('./plottingFxns')
%% Only run once after each session has been spike sorted
runMultiPostProcessing

%%
runMultiAnalysis

%%
filter = 'mec';

%%
spatialIndx = generateSpatialIndx(filter);
%%
allCells = poolAllCells(filter);

seshIndx = ismember(allCells.metadata(:,8),'ketamine');
ketamineCells = filterAllCellsStruct(allCells,seshIndx);
fprintf('done filtering ketamineCells\n');

% filter for WT mec cells with ketamine
seshIndx = ismember(ketamineCells.metadata(:,4),'WT');
fltrCells = filterAllCellsStruct(ketamineCells,seshIndx);
%%
clear allCells
clear ketamineCells
%%
%filter for WT mec cells by session
seshes = unique(cellfun(@num2str,fltrCells.metadata(:,1),'uni',0));

idxCellArray = cell(numel(seshes),1);
idxClusterArray = cell(numel(seshes),1);
decoherenceIdx = cell(numel(seshes),1);
decoherenceTimeIdx = cell(numel(seshes),1);
decoherenceTime = cell(numel(seshes),1);
decoherenceStartDelay = cell(numel(seshes),1);
Fs = [];
for i = 1:numel(seshes)
    fprintf('%i\n',i);
    seshIndx = ismember(fltrCells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(fltrCells,seshIndx);
    cells = seshCells;
    %% Get spatially binned FR over time
    all_fr_stacked = get_all_fr_stacked(cells);
    
    
    %% Calculate Downsampled FR over time
%     cellFR = cell2mat(squeeze(struct2cell(cells.FRtime)))';
%     smoothedCellFR = smoothdata(cellFR, 'sgolay',50);
% 
%     ds_factor = 100;
%     ds_sm_cellFR = downsample(smoothedCellFR, ds_factor);
%     % umap with no output
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(ds_sm_cellFR, 'verbose','none','cluster_detail','very low');
%     % umap with graphic output
%     % [reduction, umap, clusterIdentifiers, extras]=run_umap(ds_sm_cellFR, 'cluster_output','graphic','cluster_detail','very low');
%     [dch_idx, time_idx, dchTimeSec,dchStartDelaySec, Fs,trialClust] = identifyUmapDecoherenceTimeBand(clusterIdentifiers, cells, ds_factor);
%     idxClusterArray{i} = clusterIdentifiers;
%     idxCellArray{i} = trialClust;
%%     UMAP examples
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked, 'cluster_output','graphic','cluster_detail','very low');
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked, 'cluster_output','graphic','cluster_detail','very low');
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked,'python',true,'cluster_detail','very low');
    
    %% Calculate UP of spatially binned FR
    [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked, 'verbose','none','cluster_detail','very low');
    [dch_idx,dchTimeSec,dchStartDelaySec] = identifyUmapDecoherenceBand(clusterIdentifiers,cells);
    idxCellArray{i} = clusterIdentifiers;
    %%
    
    
    try
        decoherenceIdx{i} = dch_idx;
%         decoherenceTimeIdx{i} = time_idx;
    catch
        decoherenceIdx{i} = [];
    end
    
    decoherenceTime{i} = dchTimeSec;
    decoherenceStartDelay{i} = dchStartDelaySec;
end
%%
close all
% Plot Cell Array
plot_multiDimensionalCellArray(idxCellArray)
goodFigPrefs
colormap jet

% Plot Cell Array
plot_multiDimensionalCellArray(idxClusterArray)
goodFigPrefs
colormap jet

% Plot how many groups UMAP identified per session
figure()
bar(cellfun(@max,idxCellArray))
goodFigPrefs
title('Number of UMAP identified groups');

% Plot number of trials
figure();
bar(cellfun(@numel,decoherenceIdx))
goodFigPrefs
title('Decoherence Period Length (trials)');

% Highlight the Decoherence Period
highlightedDecoherenceIndx = cell(numel(seshes),1);
for i = 1:numel(seshes)
    iCA = idxCellArray{i};
    iDI = decoherenceIdx{i};
    highlightedDecoherenceIndx{i} = iCA;
    highlightedDecoherenceIndx{i}(iDI) = 10;
end
    

plot_multiDimensionalCellArray(highlightedDecoherenceIndx)
goodFigPrefs
title('Decoherence Period');
colormap jet

% Plot length of decoherence period
dchTimeMin = cell2mat(decoherenceTime)./60;
figure();
bar(dchTimeMin)
title('Decoherence Period Length (min)');
goodFigPrefs

% Plot start Delay of decoherence period
dchStartDelay = cell2mat(decoherenceStartDelay)./60;
figure();
bar(dchStartDelay)
title('Start Delay of decoherence period(min)');
goodFigPrefs

%% Plot Single Cell Raster plots with Decoherence Period Highlighted with Trial-based Dch Index
save_figs = true;
for i = 1:numel(seshes)
    close all;
    seshIndx = ismember(fltrCells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(fltrCells,seshIndx);
    plotDchRaster(decoherenceIdx,seshCells,i, save_figs)
end

%%  Plot Single Cell Raster plots with Decoherence Period Highlighted with Spatially-binned Dch Index
save_figs = false;
for i = 1:numel(seshes)
    seshIndx = ismember(fltrCells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(fltrCells,seshIndx);
    plotDchRaster_timeBinnedFR(decoherenceIdx,cells,seshID, Fs, save_figs)
end

%%
plotAllCells(allCells);