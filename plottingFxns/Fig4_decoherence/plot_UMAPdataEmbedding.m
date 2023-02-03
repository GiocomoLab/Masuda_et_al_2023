function plot_UMAPdataEmbedding(allCells,savePath,save_figs)

ds_factor = 100;
%%
seshes = unique(cellfun(@num2str,allCells.metadata(:,1),'uni',0));

close all;
plotWidth = 1500;
plotHeight = 500;
h = figure('Position',[200 200 plotWidth plotHeight]); hold on;
for i = 1:numel(seshes)
    try
        clf;
        seshIndx = ismember(allCells.metadata(:,1),seshes{i});
        cellIndx = find(seshIndx, 1, 'first');

        name = allCells.metadata{cellIndx,2};
        genotype = allCells.metadata{cellIndx,4};
        sessionDate = allCells.metadata{cellIndx,3};

        trials = downsample(allCells.trial(cellIndx).trial,ds_factor);
        reduction = allCells.dch(cellIndx).dch.umapOutput.reduction;
        clusterIdentifiers = allCells.dch(cellIndx).dch.umapOutput.clusterIdentifiers;
        if max(trials)>290
            gainStart = find(trials>290,1,'first');
        else
            gainStart = trials(end);
        end

        ax(1) = subplot(1,2,1);
        plot(downsample(reduction(:,1),5),downsample(reduction(:,2),5),'Color',[0.8,0.8,0.8],'LineWidth',0.1); hold on;
        scatter(reduction(:,1),reduction(:,2),100,trials(1:gainStart),'.');
        xlabel('UMAP X-2D')
        ylabel('UMAP Y-2D')
        goodFigPrefs
        title(sprintf('UMAP: %s%s%s%s%s Trials',name,'-',sessionDate,'-',genotype))
        colormap(ax(1),parula)
        colorbar;

        ax(2) = subplot(1,2,2);
        scatter(reduction(:,1),reduction(:,2),50,clusterIdentifiers,'.');
        xlabel('UMAP X-2D')
        ylabel('UMAP Y-2D')
        goodFigPrefs
        title(sprintf('Clusters'))
        colormap(ax(2),jet)
        set(gcf, 'Position',  [100, 100, plotWidth, plotHeight]);
        hold off;
        %%
        if save_figs
            saveas(h,fullfile(savePath,sprintf('%s%s%s%s%s.png',name,'_',sessionDate,'_',genotype)),'png');
        else
            pause
        end
    catch
       continue 
    end
end