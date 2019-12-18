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
allCellsDES = nan(count,5); % [drugFRdiff,cntrlFRdiff, drugFREffectScore, drugCorrEffectScore, spike_depth(k)]
allCellsCorrMatrix = nan(count,300,300);
allCellsCorrScoreCurve = nan(count,300);
allCellsStabilityScoreCurve = nan(count,300);
allCellsMetaData = cell(count,6);
allCellsProbeDepth = cell(count,1);
allCellsFRtime(count) = struct(); 
allCellsTrial(count) = struct(); 
allCellsSpeed(count) = struct(); 
allCellsLickT(count) = struct(); 
allCellsLickX(count) = struct(); 

sampleRate = 50; %hz
min = 15;
secInMin = 60;
minBeforeInjx = 5;
timeFRsamples = min * secInMin * sampleRate;
pre_timeFRsamples = minBeforeInjx * secInMin * sampleRate;
post_timeFRsamples = timeFRsamples - pre_timeFRsamples;

allCellsTimeFRcircaControlInjx = nan(count,timeFRsamples);
allCellsTimeFRcircaKetamineInjx = nan(count,timeFRsamples);
halfhourSampleNum = 30 * secInMin * sampleRate;
fortyfiveminSampleNum = 45 * secInMin * sampleRate;
fiftyMinSampleNum = 50 * secInMin * sampleRate;
sixtyMinSampleNum = 60 * secInMin * sampleRate;
sixtyfiveMinSampleNum = 65 * secInMin * sampleRate;
% allCellsTimeFR30minAfterKetamineInjx = nan(count,halfhourSampleNum+1);
% allCellsTimeFR45minAfterKetamineInjx = nan(count,fortyfiveminSampleNum+1);
% allCellsTimeFRneg5to45minAfterKetamineInjx = nan(count,fiftyMinSampleNum);
allCellsTimeFRneg5to60minAfterKetamineInjx = nan(count,sixtyfiveMinSampleNum);
z = 0;
fprintf('Done with the Pre-Allocation\n');
%%
for n = 1:numel(fn)
%     try
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        session_name = sessions(n).name(1:end-4);
        animalName = extractBefore(session_name,'_');
        sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
        genotype = sessionMetaData.Genotype{n};
        gender = sessionMetaData.Gender{n};
        ketamine_day = sessionMetaData.Ketamine_day(n);
        
        seshStr = sprintf('%s_%s',animalName, sessionDate);
        trackLength = 400;
        load(fullfile(matPath), 'all_fr', 'avg_all_fr', 'all_corrmatrix', 'avg_all_corrmatrix', ...
             'all_waveforms', 'cells_to_plot','spike_depth','all_drugEffectScores',...
            'trial','all_cellCorrScore','trials_corrTemplate', 'avg_all_cellCorrScore', 'avg_cell_fr',...
            'trial_ds', 'all_frTime', 'all_cellStabilityScore','speed','lickt','lickx');

        spatialIndx = ismember(cells_to_plot,allSpatialIndx.(seshStr));
       
       if ~isempty(all_fr)
           all_fr = all_fr(spatialIndx,:,:);
           all_frTime = all_frTime(spatialIndx,:,:);
           nCells = size(all_fr,1);
           allCellsMetaData(z+1:z+nCells,:,:) = repmat({session_name, animalName, sessionDate, genotype, gender, ketamine_day},nCells,1);
           if max(trial) > 300
               maxTrial = 300;
           else
               maxTrial = max(trial);
           end
           allCellsFR(z+1:z+nCells,1:maxTrial,:) = all_fr(1:nCells,1:maxTrial,1:200);
           allCellsDES(z+1:z+nCells,:,:) = all_drugEffectScores(1:nCells,1:5);
           allCellsCorrMatrix(z+1:z+nCells,1:maxTrial,1:maxTrial) = all_corrmatrix(1:nCells,1:maxTrial,1:maxTrial);
           allCellsCorrScoreCurve(z+1:z+nCells,1:maxTrial) = all_cellCorrScore(1:nCells,1:maxTrial);
           allCellsStabilityScoreCurve(z+1:z+nCells,1:maxTrial) = all_cellStabilityScore(1:nCells,1:maxTrial);
           allCellsMetaData(z+1:z+nCells,:) = repmat({session_name, animalName, sessionDate, genotype, gender, ketamine_day},nCells,1);
           
           first50indx = find(trial_ds==50, 1,'first');
           first100indx = find(trial_ds==100, 1,'first');
           allCellsTimeFRcircaControlInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first50indx-pre_timeFRsamples+1:first50indx+post_timeFRsamples);
           allCellsTimeFRcircaKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:first100indx+post_timeFRsamples);
%            allCellsTimeFR30minAfterKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx:first100indx+halfhourSampleNum);
%            if first100indx+fortyfiveminSampleNum < size(all_frTime,2)
%                allCellsTimeFR45minAfterKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx:first100indx+fortyfiveminSampleNum);
%                allCellsTimeFRneg5to45minAfterKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:first100indx+fortyfiveminSampleNum);
%            else
%                time45IndexLength = 1:size(all_frTime(1:nCells,first100indx:end),2);
%                maxTimeIndexLength = 1:size(all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:end),2);
%                
%                cellIndexLength = z+1:z+nCells;
%                allCellsTimeFR45minAfterKetamineInjx(cellIndexLength,timeIndexLength) = all_frTime(1:nCells,first100indx:end);
%                allCellsTimeFRneg5to45minAfterKetamineInjx(cellIndexLength,maxTimeIndexLength) = all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:end);
%                
%            end
           
           % Assign struct shaped data to cells
           for m = 1:nCells
               allCellsFRtime(z+m).FRtime = all_frTime(m,:); 
               allCellsTrial(z+m).trial = trial(m,:);
               allCellsSpeed(z+m).speed = speed; 
               allCellsLickT(z+m).lickt = lickt; 
               allCellsLickX(z+m).lickx = lickx;
           end
           
           
           if first100indx+sixtyMinSampleNum < size(all_frTime,2)
               allCellsTimeFRneg5to60minAfterKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:first100indx+sixtyMinSampleNum);
           else
               maxTimeIndexLength = 1:size(all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:end),2);
               cellIndexLength = z+1:z+nCells;
               allCellsTimeFRneg5to60minAfterKetamineInjx(cellIndexLength,maxTimeIndexLength) = all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:end); 
           end
            z = z + nCells;
           fprintf('Session: %d; Adding %d for %d/%d cells\n', n,nCells,z,count)
       else
           fprintf('Adding 0 for %d/%d cells\n',z,count)
       end
%     catch e
%         warning(e.message);
%         warning('FAILED: %s\n',sessions(n).name);
%     end
end
fprintf('Done with allocation\n')

%%
% SAVE ALLCELLS MAT DATA STRUCTURE
allCells.metadata = allCellsMetaData;
allCells.spatialFR = allCellsFR; 
allCells.stabilityScoreCurve = allCellsStabilityScoreCurve;
allCells.drugEffectScores = allCellsDES;
allCells.correlationMatrix = allCellsCorrMatrix;
allCells.timeFRcircaControlInjx = allCellsTimeFRcircaControlInjx;
allCells.timeFRcircaKetamineInjx = allCellsTimeFRcircaKetamineInjx;
allCells.timeFRneg5to60minAfterKetamineInjx = allCellsTimeFRneg5to60minAfterKetamineInjx;
allCells.FRtime = allCellsFRtime;
allCells.trial = allCellsTrial;
allCells.speed = allCellsSpeed;
allCells.lickT = allCellsLickT;
allCells.lickX = allCellsLickX;
%%
clearvars -except allCells