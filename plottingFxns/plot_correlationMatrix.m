function plot_correlationMatrix(allCells,cellIndx, filter)

    allCellsCorrMatrix = allCells.correlationMatrix(cellIndx,:,:);
    metadata = allCells.metadata(cellIndx,:);

    figure(); 
    imagesc(squeeze(nanmean(allCellsCorrMatrix,1)),[0 0.3])
    colorbar;
    set(gca,'FontSize',30);
    set(gca,'FontName','Helvetica');
    title(sprintf('Trial by Trial Population Correlation Matrix(%s)',filter))


    numCells = size(allCells.spatialFR,1);
    trialNum = size(allCells.spatialFR,2);
    spatialBins = size(allCells.spatialFR,3);

    flatFR2 = reshape(permute(allCells.spatialFR,[3 1 2]), [numCells*spatialBins, trialNum]);
    %P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
    corrMatrix = corr(fillmissing(flatFR2,'linear')); 
    figure(12);
    imagesc(corrMatrix); colorbar;
    set(gca,'FontSize',30);
    set(gca,'FontName','Helvetica');
    title(sprintf('Avg Cell Trial by Trial Correlation Matrix(%s)',filter))

end