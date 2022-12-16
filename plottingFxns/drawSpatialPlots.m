% sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/*.mat');
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat'); 

sessions = filterSessions(sessions, 'mec');


%%
for n = 1:numel(sessions)
    try
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        session_name = sessions(n).name(1:end-4);
        animalName = extractBefore(session_name,'_');
        sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
        trackLength = 400;
        load(fullfile(matPath), 'all_fr', 'avg_all_fr', 'all_corrmatrix', 'avg_all_corrmatrix', ...
             'all_waveforms', 'cells_to_plot','spike_depth','all_drugEffectScores',...
            'trial','all_cellCorrScore','trials_corrTemplate', 'avg_all_cellCorrScore', 'avg_cell_fr');
        
%         imgDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/img';
        imgDir = '/Users/KeiMasuda/Desktop/fkm_analysis/img';
       
        
        
       session_name = sessions(n).name(1:end-4);
       animalName = extractBefore(session_name,'_');
       sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
       seshStr = sprintf('%s_%s',animalName, sessionDate);
       spatialIndx = ismember(cells_to_plot,allSpatialIndx.(seshStr));

       all_fr = all_fr(spatialIndx,:,:);
       avg_all_fr_spatial = squeeze(nanmean(all_fr,1));
       nCells = size(all_fr,1);
        %% Plot Waveforms for Every Cell
        close all
        figure(1);
        set(gcf,'Position',[100 100 4000 1000]); 
        rows = ceil(sqrt(nCells));
        for i = 1:nCells
            subplot(rows, rows, i)
            plot(all_waveforms(i, :));
            title(sprintf('c%d, d=%d',cells_to_plot(i),round(spike_depth(i))));
        %     set(gca,'visible','off')
            set(gca,'XTick',[], 'YTick', [])
        end
        saveName = fullfile(imgDir, strcat(sessions(n).name,'_waveforms_spatial.jpg'));
        saveas(gcf,saveName, 'jpg');
        

        %% PLOT FR for each cell
        close all;
        figure(2);
        set(gcf,'Position',[100 100 1300 1000]); 
        rows = ceil(sqrt(nCells));
        %colormap('default');
        %colorbar;
        for i = 1:nCells
            subplot(rows, rows, i)
            clear g;
            x = min(trial):max(trial);
            y = squeeze(all_fr(i, :,:));
            plot(x,mean(y'));
            ttl = title(sprintf('c%d, d=%d',cells_to_plot(i),round(spike_depth(i))));
            ttl.FontSize = 10;
            set(gca,'XTick',[], 'YTick', [])
        end
        saveName = fullfile(imgDir, strcat(sessions(n).name,'_cellFR_spatial.jpg'));
        saveas(gcf,saveName, 'jpg');

        %% PLOT AVG FR FOR SESSION
        close all;
        figure(3); clf;
        set(gcf,'Position',[100 100 1300 1000]); 
        
        clear g;

        x = min(trial):max(trial);
        y = squeeze(mean(squeeze(mean(all_fr,1)),2))';

        g(1,1) = gramm('x',x,'y',y); 
        g(1,1).geom_line();
%         g(1,1) = gramm('x',x,'y',y); 
%         g(1,1).stat_summary('type','std');
        g(1,1).set_title("Ketamine's Effect on Firing Rate by Trial", 'FontSize', 30);
        g(1,1).set_names('x','Trial Number','y','Firing Rate(Hz)');
        g(1,1).set_text_options('base_size',30);

        g.draw()

        saveName = fullfile(imgDir, strcat(sessions(n).name,'_cellFR_spatial.jpg'));
        saveas(gcf,saveName, 'jpg');
        %% PLOT CORRLEATION PLOT FOR SESSION

        figure(4);
        hold on;
        p1 = plot((trials_corrTemplate:max(trial)),avg_all_cellCorrScore, 'b', 'LineWidth',0.05);
        p1.Color(4) = 0.15;
        smoothY= sgolayfilt(avg_all_cellCorrScore, 3, 25);
        plot((trials_corrTemplate:max(trial)), smoothY,'r','LineWidth',3); 
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.005 0.025]);
        set(gca,'layer','bottom');
        box off;
        axis square;
        set(gca,'FontSize',30);
        set(gca,'FontName','Helvetica');
        set(gcf,'Position',[100 100 1000 1000])
        title('Average Correlation Score')
        xlabel('Trials')
        ylabel('Correlation Score compared to First 50 Trials')

        saveName = fullfile(imgDir, strcat(sessions(n).name,'_cellCorrplot.jpg'));
        saveas(gcf,saveName, 'jpg');

    

        %% Plot Distribution of Ketamine Correlation Scores
        figure(5);
        hold on;
        histfit(all_drugEffectScores(spatialIndx,2),50, 'kernel')
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.005 0.025]);
        set(gca,'layer','bottom');
        box off;
        axis square;
        set(gca,'FontSize',30);
        set(gca,'FontName','Helvetica');
        set(gcf,'Position',[100 100 1000 1000])
        title('Distribution of Ketamine Correlation Effect Scores')
        xlabel('Pre vs Post Ketamine Correlation Effect Score')
        ylabel('Number of Cells')

        saveName = fullfile(imgDir, strcat(sessions(n).name,'_ketCorrEffectScore_spatial.jpg'));
        saveas(gcf,saveName, 'jpg');
        %% PLOT SCATTER of Corr Scorr vs Cell Depth
        figure(6);
        hold on;
        scatter(all_drugEffectScores(spatialIndx,3),all_drugEffectScores(spatialIndx,2), 150,'filled','r');
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.005 0.025]);
        set(gca,'layer','bottom');
        box off;
        axis square;
        set(gca,'FontSize',30);
        set(gca,'FontName','Helvetica');
        set(gcf,'Position',[100 100 1000 1000])
        title('Ketamine Correlation Effect Score vs Cell Depth')
        xlabel('Distance from Tip of Probe')
        ylabel('Pre vs Post Ketamine Correlation Effect Score')
        
        saveName = fullfile(imgDir, strcat(sessions(n).name,'_ketCorrEffectVScellDepth_spatial.jpg'));
        saveas(gcf,saveName, 'jpg');
        %% PLOT Scatter of Ketamine Corr Score vs FR
        figure(7);
        hold on;
        scatter(avg_cell_fr(spatialIndx),all_drugEffectScores(spatialIndx,2), 150,'filled','r');
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.005 0.025]);
        set(gca,'layer','bottom');
        box off;
        axis square;
        set(gca,'FontSize',30);
        set(gca,'FontName','Helvetica');
        set(gcf,'Position',[100 100 1000 1000])
        title('Ketamine Correlation Effect Score vs Firing Rate')
        xlabel('FR of Cell (Hz)')
        ylabel('Pre vs Post Ketamine Correlation Effect Score')
        
        saveName = fullfile(imgDir, strcat(sessions(n).name,'_ketCorrEffectVScellFR_spatial.jpg'));
        saveas(gcf,saveName, 'jpg');
        
        %% PLOT SCATTER of FR score vs Cell Depth
        figure(8);
        hold on;
        scatter(all_drugEffectScores(spatialIndx,3),all_drugEffectScores(spatialIndx,1), 150,'filled','r');
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.005 0.025]);
        set(gca,'layer','bottom');
        box off;
        axis square;
        set(gca,'FontSize',30);
        set(gca,'FontName','Helvetica');
        set(gcf,'Position',[100 100 1000 1000])
        title('Ketamine FR Effect Score vs Cell Depth')
        xlabel('Distance from Tip of Probe')
        ylabel('Ketamine FR Effect Score(postDrugFR-preDrugFR)')
        
        saveName = fullfile(imgDir, strcat(sessions(n).name,'_ketFReffectVScellDepth_spatial.jpg'));
        saveas(gcf,saveName, 'jpg');
        
        %% PLOT Scatter of Ketamine Corr Score vs ket effect on FR score
        figure(9);
        hold on;
        scatter(all_drugEffectScores(spatialIndx,1),all_drugEffectScores(spatialIndx,2), 150,'filled','r');
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.005 0.025]);
        set(gca,'layer','bottom');
        box off;
        axis square;
        set(gca,'FontSize',30);
        set(gca,'FontName','Helvetica');
        set(gcf,'Position',[100 100 1000 1000])
        title('Ketamine Correlation Effect Score vs Firing Rate Score')
        xlabel('FR Score Cell (Hz)')
        ylabel('Pre vs Post Ketamine Correlation Effect Score')
        saveName = fullfile(imgDir, strcat(sessions(n).name,'_ketFReffectVSketCorrEffect_spatial.jpg'));
        saveas(gcf,saveName, 'jpg');
        %% PLOT ket effect on FR score histogram
        figure(10);
        hold on;
        histfit(all_drugEffectScores(spatialIndx,1),50, 'kernel')
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.005 0.025]);
        set(gca,'layer','bottom');
        box off;
        axis square;
        set(gca,'FontSize',30);
        set(gca,'FontName','Helvetica');
        set(gcf,'Position',[100 100 1000 1000])
        title('Distribution of Ketamine FR Effect Scores')
        xlabel('Pre vs Post Ketamine FR Effect Score')
        ylabel('Number of Cells')

        saveName = fullfile(imgDir, strcat(sessions(n).name,'_ketFREffectScore_spatial.jpg'));
        saveas(gcf,saveName, 'jpg');
        
        %
        %% PLOT Corr Matrix for each cell
        close all;
        figure(11);
        set(gcf,'Position',[100 100 1300 1000]); 
        rows = ceil(sqrt(nCells));
        %colormap('default');
        %colorbar;
        spatialCorrMatrix = all_corrmatrix(spatialIndx,:,:);
        for i = 1:nCells
            subplot(rows, rows, i)
            imagesc(squeeze(spatialCorrMatrix(i,:,:)));
            axis off;
            box off;
        end
        saveName = fullfile(imgDir, strcat(sessions(n).name,'_cellCorrMatrix_spatial.jpg'));
        saveas(gcf,saveName, 'jpg');
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end
