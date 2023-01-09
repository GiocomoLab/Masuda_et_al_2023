function corrMatrix = calc_correlationMatrix(cells)
% Calculate
    
    numCells = size(cells.spatialFRsmooth,1);
    trialNum = size(cells.spatialFRsmooth,2);
    spatialBins = size(cells.spatialFRsmooth,3);
    smoothFR = cells.spatialFRsmooth;
    flatFR = reshape(permute(smoothFR,[3 1 2]), [numCells*spatialBins, trialNum]);
    
    %P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
    corrMatrix = corr(fillmissing(flatFR(:,1:290),'linear')); 
    
