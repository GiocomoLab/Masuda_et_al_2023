sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/*.mat');
for n = 1:numel(sessions)
    try
        [~,filename,~] = fileparts(sessions{n});
        
        [avg_all_fr, avg_all_corrmatrix] = calculateSmoothedFiringRateAllCells(matPath, trackLength);
        plot(mean(avg_all_fr,2));
        % figure(1);
        % imagesc(abs(avg_all_corrmatrix));
        % colormap('hot');
        % colorbar;
        % set(gca,'XTick',0:10:400);
        % xticklabels(xticks-100)
        % set(gca,'YTick',0:10:400);
        % yticklabels(yticks-100)
        
        % drawLicksMultiSessions(sessions{n})
        fprintf(strcat('Stitched together:', session_name,'\n'));
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions{n});
    end
end

