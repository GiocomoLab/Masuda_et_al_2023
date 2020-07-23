% run stitch sessions together, plot raster plots, combine rasters, plot
% lick data for 3 unity sessions
sessions = {... 
'/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCNj2/HCNj2_200204_keicontrasttrack_MK8011_g0',...
}; 

%'/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCNb4/HCNb4_191014_keicontrasttrack_MK8011_g0',...
% session_name = 'HCNg2_191216_keicontrasttrack_MK8011'

% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCNg2/HCNg2_191216_keicontrasttrack_MK8013_g0'
% session_name = 'HCNg2_191216_keicontrasttrack_MK8013'

lickAccuracyAllSessions = zeros(numel(sessions),300);
for n = 1:numel(sessions)
    try
       drugName = 'MK801';
       [~,filename,~] = fileparts(sessions{n});
        % session_name = 'G4_190620_keicontrasttrack_baseline+cntrlinjx+ketamine';
        % unitySessions ={'G4_190620_keicontrasttrack_baseline1', 'G4_190620_keicontrasttrack_controlinjx1','G4_190620_keicontrasttrack_ketamine1'};
        filenameParts= strsplit(filename,strcat('_',drugName));
        filename = filenameParts{1};
        session_name = strcat(filename, '_baseline+cntrlinjx+',drugName);
        
%         session_name = strcat(filename, '_baseline+',drugName);
        
        %baselineSession = dir(str cat(sessions{n},filesep,filename,'_baseline*.mat'));
        unitySessions = {strcat(filename,'_baseline1'), strcat(filename,'_controlinjx1'), strcat(filename,'_',drugName,'1')};
%         unitySessions = {strcat(filename,'_baseline1'), strcat(filename,'_',drugName,'1')};
        %Stitch file tother 
        stitchSynchedNPdata(sessions{n}, session_name, unitySessions); %unitySessions is a cell array of session names; can be a cell array of one name
        %%
        fprintf(strcat('\nStitched together:', session_name,'\n'));
        
        trackLength = 400;
        plotWidth = 160;
        plotHeight = 500;
        preDrugTrials = 100; %100
        singleSessionRasterplots(sessions{n},session_name, trackLength, plotWidth, plotHeight, preDrugTrials)
        
        size_vert = 1042; % changed from 1042 to 346
        size_horiz = 333; % changed from 333 to 667 for repeating tracks
        numrow = 6; % number of rows in final image
        combineRASTERS(session_name, size_vert, size_horiz, numrow)
        %%
        lickAccuracyAllSessions(n,:) = drawLicksSingleSessions(sessions{n},session_name);
        
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions{n});
    end
end

% close all
% %plot(mean(lickAccuracyAllSessions,1))
% close all;
% clear g;
% x = -99:200;
% y = lickAccuracyAllSessions;
% g(1,1) = gramm('x',x,'y',y);
% g(1,1).stat_summary('type','sem');
% g(1,1).set_title("Ketamine's Effect on Lick Accuracy (HCN1ko)", 'FontSize', 30);
% g(1,1).set_names('x','Trial Number','y','Lick Accuracy');
% g(1,1).set_text_options('base_size',20);
% g.draw()
