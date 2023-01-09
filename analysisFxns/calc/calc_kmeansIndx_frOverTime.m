function idx = calc_kmeansIndx_frOverTime(cells)
%% K-means clustering on a set of cell's FR over time
    
cellFR = cell2mat(squeeze(struct2cell(cells.FRtime)))';
smoothedCellFR = smoothdata(cellFR, 'sgolay',50);
%     plot(mean(smoothedCellFR'))
% imagesc(smoothedCellFR)

ds_factor = 100;
ds_sm_cellFR = downsample(smoothedCellFR, ds_factor)';

%%
%
data = ds_sm_cellFR';
eva = evalclusters(data,'kmeans','silhouette','KList',1:10);

%% K-means on Position Spikes
k_centroid = eva.OptimalK;

[idx,C] = kmeans(data,k_centroid);
imagesc(idx)


end