% re-bin rows of data so bins are [rescale] times larger
% e.g. bin_data(data, 60) to scale 1sec bins into 1min bins
function smoothed = bin_data(data, rescale)
for i = 1:size(data,1)
    smoothed(i,:) = nanmean(reshape(data(i,1:rescale*fix(size(data,2)/rescale)),rescale,[]));
end