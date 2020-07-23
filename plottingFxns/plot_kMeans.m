function plot_kMeans(cells,savefig)
nCells = size(cells.spatialFR2,1);

all_fr = cells.spatialFR10;

% trialRange = 1:max(cells.trial(1).trial);
trialRange = 1:290;
%% Get [AvgSpikesByTrial x cells]
% all_fr_trialAvg = all_fr(:,trialRange,:);
% all_fr_trialAvg = squeeze(nanmean(all_fr_trialAvg,1));
% figure(1);
% simagesc(all_fr_trialAvg);

%% Get [trials x positionFR-stackedByCells]
all_fr_stacked = all_fr(:,trialRange,:);
all_fr_stacked = permute(all_fr_stacked, [3, 1, 2]);
all_fr_stacked = reshape(all_fr_stacked,size(all_fr_stacked,1)*size(all_fr_stacked,2),[])';
% all_fr_stacked = normalize(all_fr_stacked')';
figure(1);
imagesc(all_fr_stacked);
%%
x = all_fr_stacked;
nanx = isnan(x);
t    = 1:numel(x);
x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
all_fr_stacked = x;
%%
eva = evalclusters(all_fr_stacked,'kmeans','silhouette','KList',1:10);
figure(2);
clf;
plot(eva)


%% K-means on Position Spikes
k_centroid = eva.OptimalK;
% k_centroid = 3;
[idx,C] = kmeans(all_fr_stacked,k_centroid);
figure(3); clf
scatter3(all_fr_stacked(:,1), all_fr_stacked(:,2),all_fr_stacked(:,3), 25, idx, 'filled');   % plot three clusters with different colors
colormap jet
hold on;
plot3(C(:, 1), C(:, 2), C(:, 3), 'kx');   % plot centroids

%%


%% PCA on K-means centroids
figure(4); clf
data = all_fr_stacked;
% Option 2: First reduce the dimensionality of your data using principal component analysis (PCA), and then plot the data in the principal-component space:
[standard_data, mu, sigma] = zscore(data);     % standardize data so that the mean is 0 and the variance is 1 for each variable
[coeff, score, ~]  = pca(standard_data);     % perform PCA
new_C = (C-mu)./sigma*coeff;     % apply the PCA transformation to the centroid data
scatter3(score(:, 1), score(:, 2),score(:, 3), [], idx)     % plot 3 principal components of the cluster data (three clusters are shown in different colors)
hold on
plot(new_C(:, 1), new_C(:, 2),new_C(:, 3), 'kx')     % plot 3 principal components of the centroid data
colormap jet
firstThreeComponentsAll = coeff(:, [1:3]);
projectedDataAll = data * firstThreeComponentsAll;


%% Plot Colored Spike Raster by K-means grouping
i = randi([1,nCells]);

% for i = 1:nCells
posx = cells.posX(i).posx;

spike_idx = cells.spike_idx(i);
spike_idx = spike_idx{1};
% FRtime = allCells.FRtime(i).FRtime';

trial = cells.trial(i).trial;
name = cells.metadata{i,2};
genotype = cells.metadata{i,4};
% 
figure(5);
clf;

trial_spike_idx = trial(spike_idx);
% posx_spike_idx = posx(spike_idx);
idx1 = find(idx==1);
idx2 = find(idx==2);
idx3 = find(idx==3);
group_spike_idx = ones(length(trial_spike_idx),1).*2; 
group_spike_idx(ismember(trial_spike_idx,idx1)) = 1;
group_spike_idx(ismember(trial_spike_idx,idx2)) = 2;
group_spike_idx(ismember(trial_spike_idx,idx3)) = 3; 
groupColors =  [0.2 0.2 0.2;0 0 .7;0.6,0,0]; 

x = posx(spike_idx);
y = trial(spike_idx);

trialPlotRangeIndx= y<max(trialRange);
gscatter(x(trialPlotRangeIndx),y(trialPlotRangeIndx),group_spike_idx(trialPlotRangeIndx),groupColors,'.');


% colormap('default')
set(gca, 'YDir','reverse')
ylim([0 max(trialRange)+1]);
xlim([0 400]);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000]) 
title(sprintf('Cell %d: %s,%s',i,name,genotype))
xlabel('VR cm')
ylabel('Trial Number')
legend off;
% pause
% end
%%


%% K - means on time
cellFR = cell2mat(squeeze(struct2cell(cells.FRtime)));
smoothedCellFR = smoothdata(cellFR', 'gaussian',1000);
plot(mean(smoothedCellFR'))
imagesc(smoothedCellFR)

end