function allCells = poolAllCells(filter, sessionMetaDataPath)

% sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/*.mat');
% sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat'); 
% sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_data_matrices_noSmoothing/*.mat'); 
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/combinedSesh/fr_data_matrices_noSmoothing/*.mat'); 
spatialIndexPath = '/Users/KeiMasuda/Desktop/fkm_analysis/allSpatialIndx%s.mat';

if ~exist('filter','var')
    filter = 'mec';   
end


% load(sprintf('/Users/KeiMasuda/Desktop/fkm_analysis/allSpatialIndx%s_01.mat',filter));
load(sprintf(spatialIndexPath,filter));
sessionMetaData = readtable(sessionMetaDataPath);


[sessions,sessionMetaData] = filterSessions(sessions,sessionMetaData, filter);

fn = fieldnames(allSpatialIndx);
count = 0;
seshVector = [];
for k=1:numel(fn)
    if(isnumeric(allSpatialIndx.(fn{k})))
        count = count + size(allSpatialIndx.(fn{k}),2);
        seshVector = vertcat(seshVector, repmat(k,size(allSpatialIndx.(fn{k}),2),1));
    end
end

allCellsFR2 = nan(count,300,200);
allCellsFR10 = nan(count,300,40);
allCellsDES = nan(count,5); % [drugFRdiff,cntrlFRdiff, drugFREffectScore, drugCorrEffectScore, spike_depth(k)]
allCellsCorrMatrix = nan(count,300,300);
allCellsCorrScoreCurve = nan(count,300);
allCellsStabilityScoreCurve = nan(count,300);
allCellsMetaData = cell(count,8);
allCellsProbeDepth = cell(count,1);
allCellsFRtime(count) = struct(); 
allCellsTrial(count) = struct(); 
allCellsSpeed(count) = struct();
allCellsPosT(count) = struct(); 
allCellsPosX(count) = struct(); 
allCellsLickT(count) = struct(); 
allCellsLickX(count) = struct(); 
allCellsSpikeIdx = cell(count,1);
allCellsSpatialInfo = nan(count,2);
allCellsSpatialInfoCurvesBitsPerSec = nan(count,300);
allCellsSpatialInfoCurvesBitsPerSpike = nan(count,300);
allCellsPeakiness = nan(count,300);

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

dsFactor = 50;
allCellsTimeFRneg5to60minAfterKetamineInjx = nan(count,sixtyfiveMinSampleNum/dsFactor);

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
        drug_day = sessionMetaData.drug_day(n);
        dose = sessionMetaData.dose(n);
        drug = sessionMetaData.drug{n};
        
        seshStr = sprintf('%s_%s',animalName, sessionDate);
        trackLength = 400;
        load(fullfile(matPath), 'all_fr', 'avg_all_fr', 'all_corrmatrix', 'avg_all_corrmatrix', ...
             'all_waveforms', 'cells_to_plot','spike_depth','all_drugEffectScores',...
            'trial','all_cellCorrScore','trials_corrTemplate', 'avg_all_cellCorrScore', 'avg_cell_fr',...
            'trial_ds', 'all_frTime', 'all_cellStabilityScore','post','posx','speed','lickt','lickx',...
            'all_spike_idx','all_fr10','all_spatialInfo','all_spatialInfoCurves', 'all_peakiness');
        try
            spatialIndx = ismember(cells_to_plot,allSpatialIndx.(seshStr));
        catch
            fprintf('Cannot find session in allSpatialIndx /n')
            break
        end
       
       if ~isempty(all_fr)
           all_fr = all_fr(spatialIndx,:,:);
           all_frTime = all_frTime(spatialIndx,:,:);
           nCells = size(all_fr,1);
           allCellsMetaData(z+1:z+nCells,:,:) = repmat({session_name, animalName, sessionDate, genotype, gender, drug_day, dose, drug},nCells,1);
           if max(trial) > 300
               maxTrial = 300;
           else
               maxTrial = max(trial);
           end
           allCellsFR2(z+1:z+nCells,1:maxTrial,:) = all_fr(1:nCells,1:maxTrial,1:200);
           allCellsFR10(z+1:z+nCells,1:maxTrial,:) = all_fr10(1:nCells,1:maxTrial,1:40);
           allCellsDES(z+1:z+nCells,:,:) = all_drugEffectScores(1:nCells,1:5);
           allCellsCorrMatrix(z+1:z+nCells,1:maxTrial,1:maxTrial) = all_corrmatrix(1:nCells,1:maxTrial,1:maxTrial);
           allCellsCorrScoreCurve(z+1:z+nCells,1:maxTrial) = all_cellCorrScore(1:nCells,1:maxTrial);
           allCellsStabilityScoreCurve(z+1:z+nCells,1:maxTrial) = all_cellStabilityScore(1:nCells,1:maxTrial);
           allCellsMetaData(z+1:z+nCells,:) = repmat({session_name, animalName, sessionDate, genotype, gender, drug_day, dose, drug},nCells,1);
           allCellsSpikeIdx(z+1:z+nCells,:) = all_spike_idx;
           allCellsSpatialInfo(z+1:z+nCells,:) = all_spatialInfo;

           allCellsSpatialInfoCurvesBitsPerSec(z+1:z+nCells,1:maxTrial) = squeeze(all_spatialInfoCurves(:,1:maxTrial,1));
           allCellsSpatialInfoCurvesBitsPerSpike(z+1:z+nCells,1:maxTrial) = squeeze(all_spatialInfoCurves(:,1:maxTrial,2));
           
           all_peakinessMatrix = cell2mat(all_peakiness);
           allCellsPeakiness(z+1:z+nCells,1:maxTrial) = all_peakinessMatrix(1:nCells,1:maxTrial);
           
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
               allCellsTrial(z+m).trial = trial;
               allCellsSpeed(z+m).speed = speed;
               allCellsPosT(z+m).post = post;
               allCellsPosX(z+m).posx = posx; 
               allCellsLickT(z+m).lickt = lickt; 
               allCellsLickX(z+m).lickx = lickx;
           end
           
           
           if first100indx+sixtyMinSampleNum < size(all_frTime,2)
               allCellsTimeFRneg5to60minAfterKetamineInjx(z+1:z+nCells,:) = downsample(all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:first100indx+sixtyMinSampleNum)',dsFactor)';
           else
               maxTimeIndexLength = 1:ceil(size(all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:end),2)/dsFactor);
               cellIndexLength = z+1:z+nCells;
               allCellsTimeFRneg5to60minAfterKetamineInjx(cellIndexLength,maxTimeIndexLength) = downsample(all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:end)',dsFactor)'; 
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
allCells.spatialFR2 = allCellsFR2; 
allCells.spatialFR10 = allCellsFR10;
allCells.stabilityScoreCurve = allCellsStabilityScoreCurve;
allCells.drugEffectScores = allCellsDES;
allCells.correlationMatrix = allCellsCorrMatrix;
allCells.timeFRcircaControlInjx = allCellsTimeFRcircaControlInjx;
allCells.timeFRcircaKetamineInjx = allCellsTimeFRcircaKetamineInjx;
allCells.timeFRneg5to60minAfterKetamineInjx = allCellsTimeFRneg5to60minAfterKetamineInjx;
allCells.FRtime = allCellsFRtime;
allCells.trial = allCellsTrial;
allCells.speed = allCellsSpeed;
allCells.posT = allCellsPosT;
allCells.posX = allCellsPosX;
allCells.lickT = allCellsLickT;
allCells.lickX = allCellsLickX;
allCells.spike_idx = allCellsSpikeIdx;
allCells.spatialInfo = allCellsSpatialInfo;
allCells.peakiness = allCellsPeakiness;
allCells.bitsPerSecCurve = allCellsSpatialInfoCurvesBitsPerSec;
allCells.bitsPerSpikeCurve = allCellsSpatialInfoCurvesBitsPerSpike;

%%
clearvars -except allCells

end