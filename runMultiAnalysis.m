function runMultiAnalysis(filter,combinedSessionsPath)
% Calculate and save out lickt,lickx,post,posx,speed, sp, ...
%             all_fr, avg_all_fr, all_corrmatrix, avg_all_corrmatrix,...
%             all_waveforms, cells_to_plot,spike_depth,...
%             all_drugEffectScores, trial,all_cellCorrScore,...
%             trials_corrTemplate, avg_all_cellCorrScore, avg_cell_fr,...
%             trial_ds, all_frTime,all_cellStabilityScore,all_spike_idx, all_fr10,...
%             all_spatialInfo,all_spatialInfoCurves, all_peakiness
addpath(genpath('/Users/KeiMasuda/Desktop/MalcolmFxn'));
% sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/*.mat');

% combinedSessionsPath = '/Users/KeiMasuda/Desktop/fkm_analysis/combinedSesh/*.mat';
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/combinedSesh/*.mat');
% filter = 'mec';     
sessions = filterSessions(sessions, [],filter);
%%
for n = 1:numel(sessions)
    try
        fprintf('Starting Session %i/%i...\n',n,numel(sessions));
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        trackLength = 400;
        
        [lickt,lickx,post,posx,speed, sp, ...
            all_fr, avg_all_fr, all_corrmatrix, avg_all_corrmatrix,...
            all_waveforms, cells_to_plot,spike_depth,...
            all_drugEffectScores, trial,all_cellCorrScore,...
            trials_corrTemplate, avg_all_cellCorrScore, avg_cell_fr,...
            trial_ds, all_frTime,all_cellStabilityScore,all_spike_idx, all_fr10,...
            all_spatialInfo,all_spatialInfoCurves, all_peakiness]...
            = calcFRmapCorrMatrixAllCells(matPath, trackLength);
                
%         doPCA(matPath); 

%         
        %save avg_all_fr and avg_all_corrmatrix to OAK
%         saveDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices_noSpeedFilter';
%         saveDir = '/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter';
        saveDir = '/Users/KeiMasuda/Desktop/fkm_analysis/combinedSesh/fr_data_matrices_noSmoothing';
        
        [~,sessionName,~] = fileparts(matPath);
        saveName = fullfile(saveDir, strcat(sessionName,'_fr+corr.mat'));
        save(saveName, 'all_fr', 'avg_all_fr', 'all_corrmatrix', 'avg_all_corrmatrix', ...
             'all_waveforms', 'cells_to_plot','spike_depth','all_drugEffectScores',...
            'trial','all_cellCorrScore','trials_corrTemplate', 'avg_all_cellCorrScore', 'avg_cell_fr',...
            'trial_ds', 'all_frTime','all_cellStabilityScore','post','posx','speed','lickt','lickx',...
            'all_spike_idx','all_fr10','all_spatialInfo','all_spatialInfoCurves', 'all_peakiness');     

        
        fprintf(strcat('Analyzed:', sessions(n).name,'\n'));
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end

end

