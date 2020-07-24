%%
% Only run once after each session has been spike sorted
runMultiPostProcessing

%%
% Reads sessions 
filter = 'mec';
%%
runMultiAnalysis
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
%filter for WT mec cells by session
seshes = unique(cellfun(@num2str,fltrCells.metadata(:,1),'uni',0));

idxCellArray = cell(numel(seshes),1);
decoherenceIdx = cell(numel(seshes),1);
decoherenceTime = cell(numel(seshes),1);
decoherenceStartDelay = cell(numel(seshes),1);
for i = 1:numel(seshes)
    seshIndx = ismember(fltrCells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(fltrCells,seshIndx);
    cells = seshCells;
    
    all_fr_stacked = get_all_fr_stacked(cells);
    fprintf('%i\n',i);
    close all;
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked, 'cluster_output','graphic','cluster_detail','very low');
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked, 'cluster_output','graphic','cluster_detail','very low');

%     [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked,'python',true,'cluster_detail','very low');
    [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked, 'verbose','none','cluster_detail','very low');
    dch_idx = identifyUmapDecoherenceBand(clusterIdentifiers);
    dchTimeSec = findDecoherencePeriodLength(dch_idx, cells);
    dchStartDelaySec = findDchPeriodStartDelay(dch_idx, cells);
     
%     idx = calc_kmeansIndx_frOverTime(cells);
%     idx = calc_kmeansIndx_spatialFR(cells)
    
    idxCellArray{i} = clusterIdentifiers;
    try
        decoherenceIdx{i} = dch_idx;
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

%% Plot Single Cell Raster plots with Decoherence Period Highlighted
for i = 1:numel(seshes)
    seshIndx = ismember(fltrCells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(fltrCells,seshIndx);
    plotDchRaster(decoherenceIdx,seshCells,i, save_figs)
end


%%
plotAllCells(allCells);