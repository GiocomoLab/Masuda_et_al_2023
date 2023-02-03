function plotAllSingleCells(allCells, image_save_dir, save_figs)
    addpath(genpath('./plottingFxns'))
    
    close all;
    plotWidth = 1000;
    plotHeight = 600;
    h = figure('Position',[200 200 plotWidth plotHeight]);
    for i = 1:size(allCells.spatialFR4,1)
        singleCellFR2cm = squeeze(allCells.spatialFR2(i,:,:));
        singleCellFR4cm = squeeze(allCells.spatialFR4(i,:,:));
        singleCellsmoothFR = squeeze(allCells.spatialFRsmooth(i,:,:));
        
        name = allCells.metadata{i,2};
        genotype = allCells.metadata{i,4};
        sessionDate = allCells.metadata{i,3};
%         % 
%         ogSpatialBinSize = 2;
%         spatialBinSize = 4;
%         numCol2AvgOver = spatialBinSize/ogSpatialBinSize;
%         singleCellFR = reshape(nanmean(reshape(singleCellFR4cm.',numCol2AvgOver,[])),size(singleCellFR4cm,2)/numCol2AvgOver,[]).';
% 
        posx = allCells.posX(i).posx;
        spike_idx = allCells.spike_idx(i);
        spike_idx = spike_idx{1};
        % FRtime = allCells.FRtime(i).FRtime';

        trial = allCells.trial(i).trial;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        clf;

        row = 1;
        col = 2;

        columnSubplot = [1];
        subplot(row,col,columnSubplot)
        scatter(posx(spike_idx),trial(spike_idx),'k.');
        colormap('default')
        set(gca, 'YDir','reverse')
        ylim([0 max(trial)+1]);
        xlim([0 400]);
        set(gca,'TickDir','out');
%         set(gca,'Xticklabel',[], 'Yticklabel', [])
        set(gca,'ticklength',[0.01 0.025]);   
        set(gca,'layer','bottom');
        
        set(gca,'FontName','Helvetica');
        box off;
        set(gca,'FontSize',20);
        set(gca,'FontName','Helvetica');
        % set(gcf,'Position',[100 100 1000 1000])
        title(sprintf('Cell %d: %s,%s,%s',i,name,genotype,sessionDate))
        xlabel('VR cm')
        ylabel('Trial Number')

        subplot(row,col,columnSubplot+1)
        imagesc(singleCellsmoothFR);
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.02 0.01]); 
        set(gca,'layer','bottom');
        box off;
        set(gca,'Xticklabel',[], 'Yticklabel', [])
        title(sprintf('Max FR %.1f Hz',max(singleCellsmoothFR,[],'all')))
        
   
%         set(gca,'FontSize',20);
%         set(gca,'FontName','Helvetica');
%         % set(gcf,'Position',[100 100 1000 1000])
% %         title(sprintf('Firing Rate; cellID: %d',sortIndx(i)))
%         xlabel('VR cm')
%         ylabel('Trial Number')
        %%
        if save_figs
            saveas(h,fullfile(image_save_dir,sprintf('%s%s%s%s%d.png',name,'_',sessionDate,'_',i)),'png');
        else
            pause
        end
    end
end