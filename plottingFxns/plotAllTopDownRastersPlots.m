% plot raster plots for all sessions with session at the top


sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/*/*_g0');
% Get subset of desired sessions
sessions = filterSessions(sessions, 'KO');

%%
for n = 1:numel(sessions)
    try
        drugName = 'ketamine';
        filename = sessions(n).name;
        data_dir = fullfile(sessions(n).folder,sessions(n).name);
        rmvstr = split(filename,drugName);
        session_name = strcat(rmvstr{1}, 'baseline+cntrlinjx+',drugName);
  
        trackLength = 400;
        plotWidth = 160;
        plotHeight = 500;
        preDrugTrials = 100; %100
        singleSessionRasterplots_topDown(data_dir,session_name, trackLength, plotWidth, plotHeight)
        
        size_vert = 1042; % changed from 1042 to 346
        size_horiz = 333; % changed from 333 to 667 for repeating tracks
        numrow = 6; % number of rows in final image
        combineTopDownRASTERS(session_name, size_vert, size_horiz, numrow)
        
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end