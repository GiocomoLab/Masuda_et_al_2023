function idx = calc_kmeansIndx_spatialFR(cells)
nCells = size(cells.spatialFR2,1);
all_fr = cells.spatialFR10;

% trialRange = 1:max(cells.trial(1).trial);
trialRange = 1:290;

%% Get [trials x positionFR-stackedByCells]
all_fr_stacked = all_fr(:,trialRange,:);
all_fr_stacked = permute(all_fr_stacked, [3, 1, 2]);
all_fr_stacked = reshape(all_fr_stacked,size(all_fr_stacked,1)*size(all_fr_stacked,2),[])';
all_fr_stacked = nanmean(all_fr_stacked,2);
%%
x = all_fr_stacked;
nanx = isnan(x);
t    = 1:numel(x);
x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
all_fr_stacked = x;
%%
% eva = evalclusters(all_fr_stacked,'kmeans','silhouette','KList',1:10);

%% K-means on Position Spikes
% k_centroid = eva.OptimalK;
k_centroid = 3;
[idx,~] = kmeans(all_fr_stacked,k_centroid);
end