% sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/*.mat');
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/*.mat');
sessions = filterSessions(sessions, 'mec');
%%
for n = 1:numel(sessions)
    try
        matPath = fullfile(sessions(n).folder, sessions(n).name);     
        doPCA(matPath); 
  
        fprintf(strcat('Analyzed:', sessions(n).name,'\n'));
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end

