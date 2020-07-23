function plotAllSingleCellRateMaps(allCells,save_figs)
    addpath(genpath('./plottingFxns'))
    image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/ratemaps';
%     save_figs = true;
    close all;
    plotWidth = 160;
    plotHeight = 500;
    h = figure('Position',[100 100 plotWidth plotHeight]);
    for i = 1:size(allCells.spatialFR10,1)
        singleCellFR2cm = squeeze(allCells.spatialFR2(i,:,:));

        name = allCells.metadata{i,2};
        genotype = allCells.metadata{i,4};
        sessionDate = allCells.metadata{i,3};
        % 
        ogSpatialBinSize = 2;
        spatialBinSize = 10;
        numCol2AvgOver = spatialBinSize/ogSpatialBinSize;
        singleCellFR = reshape(nanmean(reshape(singleCellFR2cm.',numCol2AvgOver,[])),size(singleCellFR2cm,2)/numCol2AvgOver,[]).';

        posx = allCells.posX(i).posx;
        spike_idx = allCells.spike_idx(i);
        spike_idx = spike_idx{1};
        % FRtime = allCells.FRtime(i).FRtime';

        trial = allCells.trial(i).trial;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clf;
% 
%         row = 4;
%         col = 2;
% 
%         columnSubplot = [1,3,5,7];
% %         subplot(row,col,columnSubplot)
%         scatter(posx(spike_idx),trial(spike_idx),'k.');
%         colormap('default')
%         set(gca, 'YDir','reverse')
%         ylim([0 max(trial)+1]);
%         xlim([0 400]);
%         set(gca,'TickDir','out');
%         set(gca,'Xticklabel',[], 'Yticklabel', [])
%         set(gca,'ticklength',[0.01 0.025]);   
%         set(gca,'layer','bottom');
%         
%         set(gca,'FontName','Helvetica');
%         box off;
% %         set(gca,'FontSize',20);
%         set(gca,'FontName','Helvetica');
%         % set(gcf,'Position',[100 100 1000 1000])
%         title(sprintf('Cell %d: %s,%s',i,name,genotype))
%         xlabel('VR cm')
%         ylabel('Trial Number')

%         subplot(row,col,columnSubplot+1)
        imagesc(singleCellFR);
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.005 0.025]);
        set(gca,'layer','bottom');
        box off;
        axis off;
%         set(gca,'FontSize',20);
        set(gca,'FontName','Helvetica');
        % set(gcf,'Position',[100 100 1000 1000])
        title(sprintf('Cell %d: %s,%s',i,name,genotype))
%         xlabel('VR cm')
%         ylabel('Trial Number')

        if save_figs
            saveas(h,fullfile(image_save_dir,sprintf('%s%s%s%s%d.png',name,'_',sessionDate,'_',i)),'png');
        else
            pause
        end
    end
end