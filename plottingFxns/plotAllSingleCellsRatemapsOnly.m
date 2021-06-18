function plotAllSingleCellsRatemapsOnly(cells, image_save_dir, save_figs)
    addpath(genpath('./plottingFxns'))
    
    close all;
    plotWidth = 160;
    plotHeight = 500;
    h = figure('Position',[200 200 plotWidth plotHeight]);
    for i = 1:size(cells.spatialFRsmooth,1)      

        singleCellsmoothFR = squeeze(cells.spatialFRsmooth(i,:,:));       
        name = cells.metadata{i,2};
        genotype = cells.metadata{i,4};
        sessionDate = cells.metadata{i,3};


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clf;
        imagesc(singleCellsmoothFR);
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.02 0.01]); 
        set(gca,'layer','bottom');
        box off;
        set(gca,'Xticklabel',[], 'Yticklabel', [])
        title(sprintf('Cell %d: %s,%s',i,name,genotype))
        colormap('hot');
%         xlabel('Virtual cm')
%         ylabel('Trial')
%         x_tick_label = get(gca,'xticklabels');
%         new_x_tick_label = cellfun(@(x) str2num(x)*2,x_tick_label);
%         xticklabels(new_x_tick_label)
        %%
        if save_figs
            saveas(h,fullfile(image_save_dir,sprintf('%s%s%s%s%d.png',name,'_',sessionDate,'_',i)),'png');
        else
            pause
        end
    end
end