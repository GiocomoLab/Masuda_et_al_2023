function dch = calcDecoherenceBandForSessions(fltrCells,dchFolderPath,paramsPath)

%filter for WT mec cells by session
seshes = unique(cellfun(@num2str,fltrCells.metadata(:,1),'uni',0));
numSesh = numel(seshes);

% load params tabel
params = readtable(paramsPath);
ds_factor = params.ds_factor;

%filter for WT mec cells by session
seshes = unique(cellfun(@num2str,fltrCells.metadata(:,1),'uni',0));

idxCellArray = cell(numel(seshes),1);
idxClusterArray = cell(numel(seshes),1);
decoherenceIdx = cell(numel(seshes),1);
decoherenceTimeIdx = cell(numel(seshes),1);
decoherenceTime = cell(numel(seshes),1);
decoherenceStartDelay = cell(numel(seshes),1);
umapArray = cell(numel(seshes),1);
Fs = [];
%%
for i = 1:numSesh
    fprintf('\nClustering session %i/%i\n',i,numel(seshes));
    seshIndx = ismember(fltrCells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(fltrCells,seshIndx);
    cells = seshCells;

    %% Calculate Downsampled FR over time
    cellFR = cell2mat(squeeze(struct2cell(cells.FRtime)))';
    smoothedCellFR = smoothdata(cellFR, 'sgolay',50);

    ds_sm_cellFR = downsample(smoothedCellFR, ds_factor);
    norm_ds_sm_cellFR = normalize(ds_sm_cellFR,2,'range');
    
    trials = downsample(cells.trial(1).trial,ds_factor);
    try
        gainTrialsStart = find(trials>290,1,'first');
        norm_ds_sm_cellFR = norm_ds_sm_cellFR(1:gainTrialsStart,:);
    catch
       fprintf('No gain trials to chop off. Less than 290 trials.\n')
    end
    
%     imagesc(norm_ds_sm_cellFR); colorbar;
    %% umap with no output
    savePath = fullfile(dchFolderPath,sprintf('umap_template_sesh%i.mat',i));
    [reduction, umap, clusterIdentifiers, extras] = run_umap(norm_ds_sm_cellFR, 'verbose','none','cluster_detail','adaptive','save_template_file', savePath); %'python',true
    
    %% umap with graphic output
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(norm_ds_sm_cellFR, 'cluster_output','graphic','n_components',2,'cluster_detail','adaptive');
    %% umap with python
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(ds_sm_cellFR, 'cluster_output','graphic','n_components',2,'cluster_detail','adaptive','python', true);
    %%
    [dch_idx, time_idx, dchTimeSec,dchStartDelaySec,Fs, trialClust,trialBasedClustIdentifiers] = identifyUmapDecoherenceTimeBand(clusterIdentifiers, cells, ds_factor);
    %%
    idxClusterArray{i} = trialBasedClustIdentifiers;
    idxCellArray{i} = trialClust;
    %%
    try
        decoherenceIdx{i} = dch_idx;
        decoherenceTimeIdx{i} = time_idx;
    catch
        decoherenceIdx{i} = [];
    end
    
    decoherenceTime{i} = dchTimeSec;
    decoherenceStartDelay{i} = dchStartDelaySec;
    umapOutput.umap = umap;
    umapOutput.reduction = reduction;
    umapOutput.clusterIdentifiers = clusterIdentifiers;
    umapOutput.extras = extras;
    umapArray{i} = umapOutput;
end  

% pass values into output struct
dch.idxCellArray = idxCellArray;
dch.idxClusterArray = idxClusterArray;
dch.decoherenceIdx = decoherenceIdx;
dch.seshes = seshes;
dch.decoherenceTime = decoherenceTime;
dch.decoherenceStartDelay = decoherenceStartDelay;
dch.decoherenceTimeIdx = decoherenceTimeIdx;
dch.Fs = Fs;
dch.umapOutput = umapArray;
dch.timeDownSample = ds_factor;
end

%% Get spatially binned FR over time
%     all_fr_stacked = get_all_fr_stacked(cells);

%% Calculate UP of spatially binned FR
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked, 'verbose','none','cluster_detail','very low');
%     [dch_idx,dchTimeSec,dchStartDelaySec] = identifyUmapDecoherenceBand(clusterIdentifiers,cells);
%     idxCellArray{i} = clusterIdentifiers;

%%     UMAP examples
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked, 'cluster_output','graphic','cluster_detail','very low');
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked, 'cluster_output','graphic','cluster_detail','very low');
%     [reduction, umap, clusterIdentifiers, extras]=run_umap(all_fr_stacked,'python',true,'cluster_detail','very low');