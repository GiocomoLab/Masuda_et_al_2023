function all_fr_stacked = get_all_fr_stacked(cells)
% Get [trials x positionFR-stackedByCells]
% Input: cells struct

nCells = size(cells.spatialFR2,1);

all_fr = cells.spatialFR10;

% trialRange = 1:max(cells.trial(1).trial);
trialRange = 1:290;
%% Get [AvgSpikesByTrial x cells]
% all_fr_trialAvg = all_fr(:,trialRange,:);
% all_fr_trialAvg = squeeze(nanmean(all_fr_trialAvg,1));
% figure(1);
% imagesc(all_fr_trialAvg);

%% Get [trials x positionFR-stackedByCells]
all_fr_stacked = all_fr(:,trialRange,:);
all_fr_stacked = permute(all_fr_stacked, [3, 1, 2]);
all_fr_stacked = reshape(all_fr_stacked,size(all_fr_stacked,1)*size(all_fr_stacked,2),[])';

end