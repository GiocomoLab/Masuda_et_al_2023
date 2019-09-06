sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/*.mat');
%%
for n = 1:numel(sessions)
    try
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        trackLength = 400;
        
        [post,posx,sp, all_fr, avg_all_fr, all_corrmatrix, avg_all_corrmatrix,...
            all_waveforms, cells_to_plot,spike_depth,...
            all_drugEffectScores, trial,all_cellCorrScore,trials_corrTemplate, avg_all_cellCorrScore, avg_cell_fr]...
            = calcFRmapCorrMatrixAllCells(matPath, trackLength);
        
        
%         doPCA(matPath); 


%         
        %save avg_all_fr and avg_all_corrmatrix to OAK
        saveDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices_noSpeedFilter';
        [~,sessionName,~] = fileparts(matPath);
        saveName = fullfile(saveDir, strcat(sessionName,'_fr+corr.mat'));
        save(saveName, 'all_fr', 'avg_all_fr', 'all_corrmatrix', 'avg_all_corrmatrix', ...
             'all_waveforms', 'cells_to_plot','spike_depth','all_drugEffectScores',...
            'trial','all_cellCorrScore','trials_corrTemplate', 'avg_all_cellCorrScore', 'avg_cell_fr');     

        
        fprintf(strcat('Analyzed:', sessions(n).name,'\n'));
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end

