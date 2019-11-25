% plot lick rasters for all sessions


sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat'); 
filter = 'WT';     
sessions = filterSessions(sessions, filter);

lickAccuracyAllSessions = zeros(numel(sessions),300);
for n = 1:numel(sessions)
    try
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        dataPath = fullfile(sessions(n).folder(1:end-30), strcat(sessions(n).name(1:end-12),'.mat'));
   

        lickAccuracyAllSessions(n,:) = drawLicksSingleSessions(sessions(n).folder(1:end-30),strcat(sessions(n).name(1:end-12)));
        fprintf(strcat('processed together:', session_name,'\n'));
        
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end

goodDuringControlIndx = size(lickAccuracyAllSessions,1);
for i = 1:size(lickAccuracyAllSessions,1)
    if mean(lickAccuracyAllSessions(i,50:100))>0.8
        goodDuringControlIndx(i) = 1;
    else
        goodDuringControlIndx(i) = 0;
    end
end

goodDuringControlIndxWT = goodDuringControlIndx;
lickAccuracyAllSessionsWT = lickAccuracyAllSessions;

%%


sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat'); 
filter = 'KO';     
sessions = filterSessions(sessions, filter);

lickAccuracyAllSessions = zeros(numel(sessions),300);
for n = 1:numel(sessions)
    try
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        dataPath = fullfile(sessions(n).folder(1:end-30), strcat(sessions(n).name(1:end-12),'.mat'));
   

        lickAccuracyAllSessions(n,:) = drawLicksSingleSessions(sessions(n).folder(1:end-30),strcat(sessions(n).name(1:end-12)));
        fprintf(strcat('processed together:', session_name,'\n'));
        
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end

goodDuringControlIndx = size(lickAccuracyAllSessions,1);
for i = 1:size(lickAccuracyAllSessions,1)
    if mean(lickAccuracyAllSessions(i,50:100))>0.8
        goodDuringControlIndx(i) = 1;
    else
        goodDuringControlIndx(i) = 0;
    end
end

goodDuringControlIndxKO = goodDuringControlIndx;
lickAccuracyAllSessionsKO = lickAccuracyAllSessions;

genotypeIndex = horzcat(repmat(0,1,size(lickAccuracyAllSessionsWT,1)), repmat(1,1,size(lickAccuracyAllSessionsKO,1)));
lickAccuracyAllSessionsCombo = vertcat(lickAccuracyAllSessionsWT,lickAccuracyAllSessionsKO);
%%
%plot(mean(lickAccuracyAllSessions,1))
close all;
clear g;
x = 1:300;
y = lickAccuracyAllSessionsCombo;
g(1,1) = gramm('x',x,'y',y,'color',genotypeIndex);
g(1,1).facet_grid([],genotypeIndex);
g(1,1).stat_summary('type','sem','setylim','true');
g(1,1).set_title("Ketamine's Effect on Lick Accuracy", 'FontSize', 30);
g(1,1).set_names('x','Trial Number','y','Lick Accuracy');
g(1,1).set_text_options('base_size',20);
g.draw()
fprintf('done')