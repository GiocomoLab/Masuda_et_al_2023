% Plot Firing Rates of Cells

addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));

%%
% sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/*.mat');
% sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices/*.mat'); 
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat'); 
imgDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/img';
% Remove strange sessions
filter = 'KO';
sessions = filterSessions(sessions, filter);

%%
trialNum = 300;
fr_mean_allCells= nan(numel(sessions),trialNum);
fr_mean_allSpatialBins = nan(numel(sessions),200);
corrmatrix_allCells = nan(numel(sessions), trialNum, trialNum);
% corrblock_allCells = nan(numel(sessions), 12, 12);

for n = 1:numel(sessions)
    try
       
       load(fullfile(sessions(n).folder, sessions(n).name));
       fr_mean_allCells(n,:) = mean(avg_all_fr(1:trialNum,:),2);
       fr_mean_allSpatialBins(n,:) = mean(avg_all_fr(1:trialNum,:),1);
       corrmatrix_allCells(n,:,:) = avg_all_corrmatrix(1:trialNum,1:trialNum);
%        corrblock_allCells(n,:,:) = avg_all_corrblock;

        session_name = sessions(n).name(1:end-4);
        animalName = extractBefore(session_name,'_');
        sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
        close all;
        clear g;
        x = -99:200;
        y = avg_all_fr';
% 
%         g(1,1) = gramm('x',x,'y',y); 
%         g(1,1).stat_summary('type','ci');
%         % g(1,1).geom_line();
%         g(1,1).axe_property('YLim',[0 30]);
%         g(1,1).set_title(sprintf("FR by Trial: %s-%s",animalName,sessionDate), 'FontSize', 40);
%         g(1,1).set_names('x','Trial Number','y','Firing Rate(Hz)');
%         g(1,1).set_text_options('base_size',20);
%         g.draw()
%         pause
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end
%%
%sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/HCN1*.mat');
% figure(1); clf;
% close all;
% clear g;
% 
% x = -99:200;
% y = fr_mean_allCells(:,1:300);
% 
% g(1,1) = gramm('x',x,'y',y); 
% g(1,1).stat_summary('type','std');
% g(1,1).geom_line();
% g(1,1).axe_property('YLim',[0 15]);
% g(1,1).set_title("Ketamine's Effect on Firing Rate by Trial", 'FontSize', 30);
% g(1,1).set_names('x','Trial Number','y','Firing Rate(Hz)');
% g(1,1).set_text_options('base_size',30);
% 
% g.draw()


%% Plot by Mouse
figure(1);clf;hold on;
if strcmp(filter, 'WT')
    mice={'G1','G2','G3','G4','G5','HCNd1','HCNe2'};
elseif strcmp(filter, 'KO')
    mice={'HCN1','HCNd2','HCNe1','HCNe3','HCNe4'};
else
    fprintf('Bad filter key. Choose: mec, WT, KO');
end


for z= 1:numel(mice)
    idx = ~cellfun('isempty',strfind({sessions.name},mice{z}));
    Y = mean(fr_mean_allCells(idx,1:300));
    interpY = fillmissing(Y,'spline');
    smoothY= sgolayfilt(interpY, 3, 25);
    mouseLine = plot(linspace(-99,200,300),smoothY,'DisplayName',sprintf('Mouse %d',z),'LineWidth',1);
%     mouseLine = plot(linspace(-99,200,300),smoothY,'DisplayName',mice{z},'LineWidth',1);
end

Y = mean(fr_mean_allCells(:,1:300));
interpY = fillmissing(Y,'spline');
smoothY= sgolayfilt(interpY, 3, 25);
plot(linspace(-99,200,300),smoothY,'k','DisplayName','Avg','LineWidth',10);


legend;
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title('Average Firing Rate by Mouse')
xlabel('Trial')
ylabel('Average FR(Hz)')


%%
figure(2);
for z= 1:numel(mice)
    clf;
    idx = ~cellfun('isempty',strfind({sessions.name},mice{z}));
    clim = [-0.1,0.2];
    imagesc(squeeze(mean(corrmatrix_allCells(idx,1:300,1:300), 1, 'omitnan')),clim);
    colormap(parula(100));
    colorbar;
    set(gca,'XTick',0:10:300);
    xticklabels(xticks-100)
    set(gca,'YTick',0:10:300);
    yticklabels(yticks-100)
    title(mice{z})
    
    saveName = fullfile(imgDir, strcat(sessions(n).name,'_avg_all_corrmatrix.jpg'));
    saveas(gcf,saveName, 'jpg');
        
end


% %% Plot By Day
% figure(5);clf;hold on;
% plot(linspace(-99,200,300),mean(fr_mean_allCells([1,10,15,18],1:300)),'DisplayName','Day 1')
% plot(linspace(-99,200,300),mean(fr_mean_allCells([2,5,11],1:300)),'DisplayName','Day 2')
% plot(linspace(-99,200,300),mean(fr_mean_allCells([3,6,12,16],1:300)),'DisplayName','Day 3')
% plot(linspace(-99,200,300),mean(fr_mean_allCells([4,7,13,17],1:300)),'DisplayName','Day 4')
% plot(linspace(-99,200,300),mean(fr_mean_allCells([8,14],1:300)),'DisplayName','Day 5')
% legend;
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('Average Firing Rate by Day')
% xlabel('Trial')
% ylabel('Average FR(Hz)')
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);

