function k_centroid = optimize_k(cells)
all_fr = cells.spatialFR2;
trialRange = 1:max(cells.trial(1).trial);
% 
% %% Get [AvgSpikesByTrial x cells]
% all_fr_trialAvg = all_fr(:,trialRange,:);
% all_fr_trialAvg = squeeze(nanmean(all_fr_trialAvg,1));
% figure(1);
% imagesc(all_fr_trialAvg);

%% Get [trials x positionFR-stackedByCells]
all_fr_stacked = all_fr(:,trialRange,:);
all_fr_stacked = permute(all_fr_stacked, [3, 1, 2]);
all_fr_stacked = reshape(all_fr_stacked,size(all_fr_stacked,1)*size(all_fr_stacked,2),max(trialRange))';

eva = evalclusters(all_fr_stacked,'kmeans','CalinskiHarabasz','KList',1:10);
k_centroid = eva.OptimalK;
%% optimize with temporally binned FR
cellFR = cell2mat(squeeze(struct2cell(cells.FRtime)));
% smoothedCellFR = smoothdata(cellFR, 'gaussian',10);
% plot(smoothedCellFR(1,1:500))
eva = evalclusters(smoothedCellFR,'kmeans','CalinskiHarabasz','KList',1:10);
plot(eva)
k_centroid = eva.OptimalK;

end