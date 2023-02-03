function plot_cellClustersOverTime(cells,savefig)
%% K-means clustering on a set of cell's FR over time
    
cellFR = cell2mat(squeeze(struct2cell(cells.FRtime)))';
smoothedCellFR = smoothdata(cellFR, 'sgolay',50);
%     plot(mean(smoothedCellFR'))
imagesc(smoothedCellFR)

ds_factor = 100;
ds_sm_cellFR = downsample(smoothedCellFR, ds_factor)';
figure(1);clf;
imagesc(ds_sm_cellFR)
%%
%
data = ds_sm_cellFR';
eva = evalclusters(data,'kmeans','silhouette','KList',1:10);
figure(2);
clf;
plot(eva)


%% K-means on Position Spikes
k_centroid = eva.OptimalK;
% k_centroid = 3;
[idx,C] = kmeans(data,k_centroid);
figure(3); clf
scatter3(data(:,1), data(:,2),data(:,3), 25, idx, 'filled');   % plot three clusters with different colors
colormap jet
hold on;
plot3(C(:, 1), C(:, 2), C(:, 3), 'kx');   % plot centroids

%% PCA on K-means centroids
figure(4); clf
% Option 2: First reduce the dimensionality of your data using principal component analysis (PCA), and then plot the data in the principal-component space:
[standard_data, mu, sigma] = zscore(data);     % standardize data so that the mean is 0 and the variance is 1 for each variable
[coeff, score, ~]  = pca(standard_data);     % perform PCA
new_C = (C-mu)./sigma*coeff;     % apply the PCA transformation to the centroid data
scatter3(score(:, 1), score(:, 2),score(:, 3), [], idx)     % plot 2 principal components of the cluster data (three clusters are shown in different colors)
hold on
plot(new_C(:, 1), new_C(:, 2),new_C(:, 3), 'kx')     % plot 2 principal components of the centroid data
colormap jet
firstThreeComponentsAll = coeff(:, [1:3]);
projectedDataAll = data * firstThreeComponentsAll;

%% Plot kmeans Indx
figure(5); clf
imagesc(idx)

end