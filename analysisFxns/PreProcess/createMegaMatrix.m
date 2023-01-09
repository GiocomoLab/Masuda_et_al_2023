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
        seshVector = vertcat(seshVector, repmat(k,size(allSpatialIndx.(fn{k}),2),1));
    end
end

allCellsFR = nan(count,300,200);

metaData = struct;
for i = 1:count
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

sz = [size(allCellsFR,1)*size(allCellsFR,2), 10];
varTypes = {'string','double','double','double','string','string','double','double','double','logical'};
varNames = {'animalName', 'sessionDate', 'trialNum','cellNum','gender', ...
    'genotype','weight_g', 'ketamine_day', 'timeSinceKetamine','ketamineAdministered'};

% varNames = {'animalName', 'sessionDate', 'trialNum','cellNum','gender', ...
%     'genotype','weight_g', 'ketamine_day','trialStartTime','ketamineInjxTimeSec', 'controlInjxTimeSec'};
allTrialsAllCellsTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

%%
z = 0;
y = 0;
samplingRate = 50; %Hz
for n = 1:numel(fn)
    matPath = fullfile(sessions(n).folder, sessions(n).name);
    session_name = sessions(n).name(1:end-4);
    animalName = extractBefore(session_name,'_');
    sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
    seshStr = sprintf('%s_%s',animalName, sessionDate);
    trackLength = 400;
    load(fullfile(matPath), 'all_fr', 'cells_to_plot','spike_depth','trial');
    spatialIndx = ismember(cells_to_plot,allSpatialIndx.(seshStr));

    if ~isempty(all_fr)
        all_fr = all_fr(spatialIndx,:,:);
        nCells = size(all_fr,1);
        allCellsFR(z+1:z+nCells,1:max(trial),:) = all_fr(1:nCells,1:max(trial),1:200);

        ketamineInjxTimeSec = find(trial==100,1)/samplingRate;
        controlInjxTimeSec = find(trial==50,1)/samplingRate;

        
        
        for m = 1:nCells
            mdIndx = z+1;
            metaData(mdIndx).ketamineInjxTimeSec = ketamineInjxTimeSec;
            metaData(mdIndx).controlInjxTimeSec = controlInjxTimeSec;
            animalName = string(metaData(mdIndx).animalName);
            sessionDate = metaData(mdIndx).sessionDate;
            gender = metaData(mdIndx).gender;
            genotype = metaData(mdIndx).genotype;
            weight_g = metaData(mdIndx).weight_g;
            ketamine_day = metaData(mdIndx).ketamine_day;
            
            rows = y+1:(y + max(trial));
            trialRange = 1:max(trial);
            trialStartTime = arrayfun(@(x) find(trial==x,1)/samplingRate,trialRange','UniformOutput',false);
            timeSinceKetamine = cell2mat(trialStartTime)-ketamineInjxTimeSec;
            ketamineAdminBool = timeSinceKetamine>=0;
            
            rowCount = numel(rows);
            rowData = table(...
                repmat(animalName,rowCount,1),...
                repmat(str2double(sessionDate),rowCount,1),...
                trialRange',...
                repmat(m,rowCount,1),...
                repmat(gender,rowCount,1),...
                repmat(genotype,rowCount,1),...
                repmat(weight_g,rowCount,1),...
                repmat(ketamine_day,rowCount,1),...
                timeSinceKetamine,...
                ketamineAdminBool...
            );
            

            allTrialsAllCellsTable(rows,:) = rowData;
          
            y = y + max(trial);
%             fprintf('Trial: %d for session %d\n', y,n)
        end

        z = z + nCells;
        fprintf('Session: %d; Adding %d for %d/%d cells\n', n,nCells,z,count)
    else
       fprintf('Adding 0 for %d/%d cells\n',z,count)
    end
end

fprintf('done\n')

%%
spatialBin = reshape(allCellsFR(:,1:300,1:200),[],size(allCellsFR,3));
%%
numCol2AvgOver = 5;
spatialBin10cm = reshape(nanmean(reshape(spatialBin.',numCol2AvgOver,[])),size(spatialBin,2)/numCol2AvgOver,[]).';

allT_allC_frTable = array2table(spatialBin10cm);
% Remove Nans
allT_allC_frTable = fillmissing(allT_allC_frTable,'previous');

%%
% Cocncatentate into megatbale
megaTable = horzcat(allTrialsAllCellsTable,allT_allC_frTable);
%%
ketBoolVector = allTrialsAllCellsTable(:,10);

%%
writetable(megaTable,"/Users/KeiMasuda/Desktop/allTrialsAllCellsTable10cm.csv");



%%
writetable(allT_allC_frTable, "/Users/KeiMasuda/Desktop/neuralData10cm.csv");
writetable(ketBoolVector, "/Users/KeiMasuda/Desktop/ketBool.csv");
