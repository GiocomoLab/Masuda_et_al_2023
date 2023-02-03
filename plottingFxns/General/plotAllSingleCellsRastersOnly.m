function plotAllSingleCellsRastersOnly(cells, image_save_dir, save_figs)
    addpath(genpath('./plottingFxns'))
    
    close all;
    plotWidth = 160;
    plotHeight = 500;
    h = figure('Position',[200 200 plotWidth plotHeight]);
    for i = 1:size(cells.spatialFR4,1)

        name = cells.metadata{i,2};
        genotype = cells.metadata{i,4};
        sessionDate = cells.metadata{i,3};

        posx = cells.posX(i).posx;
        spike_idx = cells.spike_idx(i);
        spike_idx = spike_idx{1};
      
        trial = cells.trial(i).trial;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clf;

        scatter(posx(spike_idx),trial(spike_idx),'k.'); hold on;
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
        title(sprintf('C%d: %s,%s',i,name,genotype))
        
   
        %%
        if save_figs
            saveas(h,fullfile(image_save_dir,sprintf('%s%s%s%s%d.png',name,'_',sessionDate,'_',i)),'png');
        else
            pause
        end
    end
end