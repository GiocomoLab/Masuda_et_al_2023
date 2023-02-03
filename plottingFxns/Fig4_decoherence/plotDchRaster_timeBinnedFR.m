function plotDchRaster_timeBinnedFR(decoherenceIdx,cells,seshID, Fs, save_figs)

    addpath(genpath('./plottingFxns'))
    image_save_dir = '/Users/KeiMasuda/Desktop/fkm_analysis/rasters_dch_timeBinned';

%%
	dchRange = decoherenceIdx{seshID};
    
%%
    plotWidth = 160;
    plotHeight = 500;
    h = figure('Position',[100 100 plotWidth plotHeight]); hold on;
    for i = 1:size(cells.spatialFR10,1)

        name = cells.metadata{i,2};
        genotype = cells.metadata{i,4};
        sessionDate = cells.metadata{i,3};
        

        posx = cells.posX(i).posx;
        post = cells.posT(i).post;
        spike_idx = cells.spike_idx(i);
        spike_idx = spike_idx{1};
        trial = cells.trial(i).trial;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clf;

        scatter(posx(spike_idx),trial(spike_idx),'k.'); hold on;
        
        x = posx(spike_idx);
        y = trial(spike_idx);
        z = post(spike_idx);
        dchStart = min(dchRange);
        dchEnd = max(dchRange);
        dchRange_timeInx = z>dchStart & z<dchEnd;
        x = x(dchRange_timeInx);
        y = y(dchRange_timeInx);
        scatter(x,y,'r.');
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


