function spatialIndx = generateSpatialIndx(filter)


addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));

%%
if ~exist('filter','var')
    filter = 'mec';   
end
% sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/*.mat');
% sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat');
% sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_data_matrices_noSmoothing/*.mat');
sessionsPath = '/Users/KeiMasuda/Desktop/fkm_analysis/combinedSesh/fr_data_matrices_noSmoothing/*.mat';
sessions = dir(sessionsPath);
% Get subset of desired sessions
% filter = 'mec';
sessions = filterSessions(sessions, [],filter);

%%
totalCell = 0;
totalSpatialCell = 0;
allSpatialIndx = {};
for n = 1:numel(sessions)
    try
        
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        session_name = sessions(n).name(1:end-4);
        animalName = extractBefore(session_name,'_');
        sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
        trackLength = 400;
        load(fullfile(matPath), 'cells_to_plot');
        
        nCells = numel(cells_to_plot);
        totalCell = totalCell + nCells;
        
        spatialIndx = [];

        for k = 1:numel(cells_to_plot)       
            totalSpatialCell = totalSpatialCell + 1;
            spatialIndx(size(spatialIndx,2)+1) = cells_to_plot(k);
        end
        
        seshIndx = sprintf('%s_%s',animalName,sessionDate);
        allSpatialIndx.(seshIndx) = spatialIndx;

    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end
fprintf(strcat('\nTotal Number of Cells: ',num2str(totalCell),'\nNumber of Filtered Cells:',num2str(totalSpatialCell),'\n'));
%%
save(sprintf('/Users/KeiMasuda/Desktop/fkm_analysis/allSpatialIndx%s.mat',filter),'allSpatialIndx');

end