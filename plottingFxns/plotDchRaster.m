function plotDchRaster(dchRange,cells, save_figs,image_save_dir)

    addpath(genpath('./plottingFxns'))
    image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/rasters_dch';

    
%%
    plotWidth = 160;
    plotHeight = 500;
    h = figure('Position',[100 100 plotWidth plotHeight]); hold on;
    for i = 1:size(cells.spatialFR2,1)
        
        
        singleCellFR2cm = squeeze(cells.spatialFR2(i,:,:));

        name = cells.metadata{i,2};
        genotype = cells.metadata{i,4};
        sessionDate = cells.metadata{i,3};
        % 
        ogSpatialBinSize = 2;
        spatialBinSize = 10;
        numCol2AvgOver = spatialBinSize/ogSpatialBinSize;
        singleCellFR = reshape(nanmean(reshape(singleCellFR2cm.',numCol2AvgOver,[])),size(singleCellFR2cm,2)/numCol2AvgOver,[]).';

        posx = cells.posX(i).posx;
        spike_idx = cells.spike_idx(i);
        spike_idx = spike_idx{1};
        % FRtime = allCells.FRtime(i).FRtime';

        trial = cells.trial(i).trial;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clf;

        scatter(posx(spike_idx),trial(spike_idx),'k.'); hold on;
        
        if ~isempty(dchRange)
            x = posx(spike_idx);
            y = trial(spike_idx);
            x = x(y>min(dchRange) & y<max(dchRange));
            y = y(y>min(dchRange) & y<max(dchRange));
            scatter(x,y,'r.');
        end
        colormap('default')
        
        set(gca, 'YDir','reverse')
        ylim([0 max(trial)+1]);
        xlim([0 400]);
        set(gca,'TickDir','out');
        set(gca,'Xticklabel',[], 'Yticklabel', [])
        set(gca,'ticklength',[0.01 0.025]);   
        set(gca,'layer','bottom');
        
        set(gca,'FontName','Helvetica');
        box off;

        set(gca,'FontName','Helvetica');
        % set(gcf,'Position',[100 100 1000 1000])
        title(sprintf('Cell %d: %s,%s',i,name,genotype))
        
        if save_figs
            saveas(h,fullfile(image_save_dir,sprintf('%s%s%s%s%d.png',name,'_',sessionDate,'_',i)),'png');
        else
            pause
            
        end

    end
end


