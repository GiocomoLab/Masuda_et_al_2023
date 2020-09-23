function plot_UMAPdataEmbedding(allCells,savePath,save_figs)


seshes = unique(cellfun(@num2str,allCells.metadata(:,1),'uni',0));

close all;
plotWidth = 600;
plotHeight = 600;
h = figure('Position',[200 200 plotWidth plotHeight]);
for i = 1:numel(seshes)
    seshIndx = ismember(allCells.metadata(:,1),seshes{i});
    cellIndx = find(seshIndx, 1, 'first');
    
    name = allCells.metadata{cellIndx,2};
    genotype = allCells.metadata{cellIndx,4};
    sessionDate = allCells.metadata{cellIndx,3};
    
    
    reduction = allCells.dch(cellIndx).dch.umapOutput.reduction;
    clusterIdentifiers = allCells.dch(cellIndx).dch.umapOutput.clusterIdentifiers;
    scatter(reduction(:,1),reduction(:,2),20,clusterIdentifiers,'.')
    colormap(jet)
    xlabel('UMAP X-2D')
    ylabel('UMAP Y-2D')
    goodFigPrefs
    title(sprintf('UMAP: %s%s%s%s%s',name,'-',sessionDate,'-',genotype))
    %%
    if save_figs
        saveas(h,fullfile(savePath,sprintf('%s%s%s%s%s.png',name,'_',sessionDate,'_',genotype)),'png');
    else
        pause
    end
end