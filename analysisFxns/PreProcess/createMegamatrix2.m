% Make allTrialsAllCellsMatrix

% sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/*.mat');
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat'); 
filter = 'mec';     
sessions = filterSessions(sessions, filter);
% load(sprintf('/Users/KeiMasuda/Desktop/fkm_analysis/allSpatialIndx%s_01.mat',filter));
load(sprintf('/Users/KeiMasuda/Desktop/fkm_analysis/allSpatialIndx%s.mat',filter));

sessionMetaData = readtable('/Users/KeiMasuda/Desktop/fkm_analysis/SessionList.xlsx');

fn = fieldnames(allSpatialIndx);
count = 0;
seshVector = [];
for k=1:numel(fn)
    if(isnumeric(allSpatialIndx.(fn{k})))
        count = count + size(allSpatialIndx.(fn{k}),2);
        seshVector = vertcat(seshVector, k);
    end
end

allCellsFR = nan(count,100,200);

metaData = struct;
for i = 1:numel(fn)
    session_name = fn{seshVector(i)};
    animalName = extractBefore(session_name,'_');
    sessionDate = extractAfter(session_name,'_');
    
    metaData(i).sessionName = session_name;
    metaData(i).animalName = animalName;
    metaData(i).sessionDate = sessionDate;
    metaData(i).gender = string(sessionMetaData(seshVector(i),:).Gender);
    metaData(i).genotype = string(sessionMetaData(seshVector(i),:).Genotype);
    metaData(i).weight_g = sessionMetaData(seshVector(i),:).Weight_g_;
    metaData(i).ketamine_day = sessionMetaData(seshVector(i),:).Ketamine_day;
end
fprintf('Done with the Pre-Allocation\n');
%%
varTypes = {'string','double','double','double','string','string','double','double',...
    'double','double','double','double','double','double','double','double','double','double','logical'};
varNames = {'animalName','sessionDate','trialNum','totalCellNum','gender','genotype','weight_g',...
    'ketamine_day','correlationScore','lickAccuracy','lickNumber','avgFR','avgSingleCellVariance',...
    'varianceFR','avgTrialSpeed','varianceSpeed','medianCellDepth','timeSinceKetamine','ketamineAdministered'};
sz = [numel(fn)*size(allCellsFR,2), numel(varTypes)];
sessionTrialTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

%%
z = 0;
y = 0;
samplingRate = 50; %Hz
% trialRange = 51:150;
trialRange = 101:200;
numTrialsForTable = numel(trialRange);
for n = 1:numel(fn)
    matPath = fullfile(sessions(n).folder, sessions(n).name);
    dataPath = fullfile(sessions(n).folder(1:end-30), strcat(sessions(n).name(1:end-12),'.mat'));
    session_name = sessions(n).name(1:end-4);
    animalName = extractBefore(session_name,'_');
    sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
    seshStr = sprintf('%s_%s',animalName, sessionDate);
    trackLength = 400;
    load(fullfile(matPath), 'all_fr', 'cells_to_plot','trial','all_cellCorrScore','spike_depth');
    load(dataPath, 'lickt','lickx', 'post','posx')
    
    spatialIndx = ismember(cells_to_plot,allSpatialIndx.(seshStr));

    all_fr = all_fr(spatialIndx,trialRange,:);
    
    nCells = size(all_fr,1);
    allCellsFR(z+1:z+nCells,1:numTrialsForTable,:) = all_fr(1:nCells,1:numTrialsForTable,1:200);

    ketamineInjxTimeSec = find(trial==100,1)/samplingRate;
    controlInjxTimeSec = find(trial==50,1)/samplingRate;

        
    animalName = string(metaData(n).animalName);
    sessionDate = metaData(n).sessionDate;
    gender = metaData(n).gender;
    genotype = metaData(n).genotype;
    weight_g = metaData(n).weight_g;
    ketamine_day = metaData(n).ketamine_day;

    rows = y+1:(y + numTrialsForTable);

    trialStartTime = arrayfun(@(x) find(trial==x,1)/samplingRate,trialRange','UniformOutput',false);
    timeSinceKetamine = cell2mat(trialStartTime)-ketamineInjxTimeSec;
    ketamineAdminBool = timeSinceKetamine>=0;


    [~,~,lick_idx] = histcounts(lickt,post);

    lickAccuracyByTrial = zeros(1,numTrialsForTable);
    lickNumByTrial = zeros(1,numTrialsForTable);

    for i = 1:numTrialsForTable
        j = trialRange(i);
        trialLicks = lickx(trial(lick_idx) == j);
        goodLicks = sum(trialLicks<5) + sum(trialLicks>max(posx)-15); 
        if trialLicks ~= 0
            lickAccuracyByTrial(i) = goodLicks/numel(trialLicks);
            lickNumByTrial(i) = numel(trialLicks);
        else
            lickAccuracyByTrial(i) = 0.0;
            lickNumByTrial(i) = 0.0;
        end
    end

    corrScoreByTrial = nanmean(all_cellCorrScore(:,trialRange),1); %avg corr scorr for each trial avg across all cells

    avgFRbyTrial = nanmean(squeeze(nanmean(all_fr,1)),2);
    avgSingleCellVariance = var(squeeze(nanmean(all_fr,3)),0,1);
    varianceFR = var(squeeze(nanmean(all_fr,1)),0,2);

    trialSpeed = arrayfun(@(x) 400/(numel(find(trial==x))/samplingRate),trialRange','UniformOutput',false);
    trialSpeed = cell2mat(trialSpeed); %speed in cm/s

    varianceSpeed = arrayfun(@(x) var(diff(posx(find(trial==x))))*samplingRate,trialRange','UniformOutput',false);
    varianceSpeed = cell2mat(varianceSpeed);

    medianCellDepth = median(spike_depth);

    rowCount = numel(rows);
    rowData = table(...
        repmat(animalName,rowCount,1),...
        repmat(str2double(sessionDate),rowCount,1),...
        trialRange',...
        repmat(numel(cells_to_plot),rowCount,1),...
        repmat(gender,rowCount,1),...
        repmat(genotype,rowCount,1),...
        repmat(weight_g,rowCount,1),...
        repmat(ketamine_day,rowCount,1),...
        corrScoreByTrial',...
        lickAccuracyByTrial',...
        lickNumByTrial',...
        avgFRbyTrial,...
        avgSingleCellVariance',...
        varianceFR,...
        trialSpeed,...
        varianceSpeed,...
        repmat(medianCellDepth,rowCount,1),...
        timeSinceKetamine,...
        ketamineAdminBool...
    );


    sessionTrialTable(rows,:) = rowData;
    y = y + numTrialsForTable;


    z = z + nCells;
    fprintf('Session: %d; Adding %d for %d/%d cells\n', n,nCells,z,count)

end

fprintf('done\n')


%%
% writetable(sessionTrialTable,"/Users/KeiMasuda/Dropbox/1_SMS/NEURS/Electives/MS&E 226/Project/sessionTrialTable.csv");
writetable(sessionTrialTable,"/Users/KeiMasuda/Dropbox/1_SMS/NEURS/Electives/MS&E 226/Project/postKetamineTable.csv");
% [row, col] = find(ismissing(sessionTrialTable))

% %%
% postKetamineTrials = sessionTrialTable(sessionTrialTable.trialNum >100,:);
% writetable(postKetamineTrials,"/Users/KeiMasuda/Dropbox/1_SMS/NEURS/Electives/MS&E 226/Project/postKetamineTable.csv");

%%
WTtable = sessionTrialTable(sessionTrialTable.genotype == "WT",:);
HCNkoTable = sessionTrialTable(sessionTrialTable.genotype == "KO",:);
%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
close all;
clear g;
g = gramm('x',categorical(sessionTrialTable.genotype),'y',sessionTrialTable.avgTrialSpeed,'color',categorical(sessionTrialTable.genotype));
% g.facet_grid([],categorical(sessionTrialTable.gender))
g.set_names('x','Trial','y','lickAccuracy', 'column','Gender')
% g.stat_summary()
g.stat_violin()
g.draw();
