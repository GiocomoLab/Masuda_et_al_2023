function plot_clustergram(cells)
    nCells = size(cells.spatialFR2,1);
    
    all_fr = [];
    ds_factor = 100;
    
    trial = squeeze(cell2mat(struct2cell(cells.trial)));
    trial = trial(:,1);
    trial_ds = downsample(trial, ds_factor); 
    
    
    
    % gaussian filter for smoothing
    smoothSigma = 10;
    smoothWindow = floor(smoothSigma*5/2)*2+1;
    gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);
    
    % Calculating smoothed firing rate
    for i = 1:nCells
        cellFR = cells.FRtime(i).FRtime;
        % smooth firing rate
        fr_smoothed = conv(repmat(cellFR,1,3),gauss_filter,'same');
        fr = fr_smoothed(numel(cellFR)+1:numel(cellFR)*2);
        %
        all_fr(i, :) = fr;
    end 
    
    all_fr = all_fr';
    all_fr_ds = downsample(all_fr, ds_factor);
    
    cgo_all=clustergram(all_fr_ds,'Standardize','Column','RowLabels',trial_ds','ColumnLabels',1:106,'Colormap',redbluecmap)
    
    plot(str2num(cell2mat(cgo_all.RowLabels)))
    colorbar

end