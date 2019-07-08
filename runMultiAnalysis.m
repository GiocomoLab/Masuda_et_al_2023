sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/*.mat');
for n = 1:numel(sessions)
    try
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        trackLength = 400;
        [avg_all_fr, avg_all_corrmatrix] = calculateSmoothedFiringRateAllCells(matPath, trackLength);
        figure(1);
        plot(mean(avg_all_fr,2));
        
        figure(2);
        imagesc(abs(avg_all_corrmatrix));
        colormap('default');
        colorbar;
        set(gca,'XTick',0:10:400);
        xticklabels(xticks-100)
        set(gca,'YTick',0:10:400);
        yticklabels(yticks-100);
        
        %save avg_all_fr and avg_all_corrmatrix to OAK
        saveDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices';
        saveName = fullfile(saveDir, strcat(sessions(n).name,'_fr+corr.mat'));
        save(saveName, 'avg_all_fr','avg_all_corrmatrix');

        % drawLicksMultiSessions(sessions{n})
        fprintf(strcat('Analyzed:', sessions(n).name,'\n'));
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions{n});
    end
end

