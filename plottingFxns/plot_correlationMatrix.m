function plot_correlationMatrix(cells,filter)

    allCellsCorrMatrix = cells.correlationMatrix;
    metadata = cells.metadata;

    figure(); 
    imagesc(squeeze(nanmean(allCellsCorrMatrix,1)),[0 0.2])
%     imagesc(squeeze(nanmean(allCellsCorrMatrix,1)))
    colorbar;
    set(gca,'FontSize',30);
    set(gca,'FontName','Helvetica');
    title(sprintf('Avg Cell Correlation Matrix(%s)',filter))
    box off;
    axis square;
    set(gca,'TickDir','out');
    set(gca,'ticklength',[0.005 0.025]);
    set(gca,'layer','bottom');  

    
%     [corrMatrix,pval] = corr(squeeze(nanmean(cells.spatialFR2,3))); 
    
    numCells = size(cells.spatialFR10,1);
    trialNum = size(cells.spatialFR10,2);
    spatialBins = size(cells.spatialFR10,3);
    flatFR2 = reshape(permute(cells.spatialFR10,[3 1 2]), [numCells*spatialBins, trialNum]);
    %P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
    corrMatrix = corr(fillmissing(flatFR2,'linear')); 
    figure();
    imagesc(corrMatrix); colorbar;
%     imagesc(corrMatrix,[0.5 0.8]); colorbar;
    set(gca,'FontSize',30);
    set(gca,'FontName','Helvetica');
    title(sprintf('Population Correlation Matrix(%s)',filter))
    box off;
    axis square;
    set(gca,'TickDir','out');
    set(gca,'ticklength',[0.005 0.025]);
    set(gca,'layer','bottom');  
end