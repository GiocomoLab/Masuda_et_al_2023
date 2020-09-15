function plotDecoherenceBandStats(dch)

close all
%% Plot Cluster Groups by Trials Array
plot_multiDimensionalCellArray(dch.idxCellArray)
goodFigPrefs
colormap jet
title('Cluster Groups by Trials Array');

%% Highlight the Decoherence Period for Trials
highlightedDecoherenceIndx = cell(numel(dch.seshes),1);
for i = 1:numel(dch.seshes)
    iCA = dch.idxCellArray{i};
    iDI = dch.decoherenceIdx{i};
    highlightedDecoherenceIndx{i} = iCA;
    highlightedDecoherenceIndx{i}(iDI) = 10;
end
plot_multiDimensionalCellArray(highlightedDecoherenceIndx)
goodFigPrefs
title('Decoherence Period');
colormap jet  

%%  Plot Cluster Groups by Time Array
plot_multiDimensionalCellArray(dch.idxClusterArray)
goodFigPrefs
colormap jet
title('Cluster Groups by Time Array');
colormap jet
yt = get(gca, 'YTick'); 
set(gca, 'YTick',yt, 'YTickLabel',round(yt*dch.Fs/60))
ylabel('Minutes');

%% Highlight the Decoherence Period for Time
highlightedDecoherenceIndx = cell(numel(dch.seshes),1);
for i = 1:numel(dch.seshes)
    iCA = dch.idxClusterArray{i};
    iDI = dch.decoherenceTimeIdx{i};
    highlightedDecoherenceIndx{i} = iCA;
    highlightedDecoherenceIndx{i}(min(iDI):max(iDI)) = 15;
end

plot_multiDimensionalCellArray(highlightedDecoherenceIndx)
goodFigPrefs
title('Decoherence Period');
colormap jet
yt = get(gca, 'YTick'); 
set(gca, 'YTick',yt, 'YTickLabel',round(yt*dch.Fs/60))
ylabel('Minutes');
%% Plot how many groups UMAP identified per session
figure()
bar(cellfun(@max,dch.idxCellArray))
goodFigPrefs
title('Number of UMAP identified groups');

%% Plot number of trials
figure();
bar(cellfun(@numel,dch.decoherenceIdx))
goodFigPrefs
title('Decoherence Period Length (trials)');

%% Plot length of decoherence period
dchTimeMin = cell2mat(dch.decoherenceTime)./60;
figure();
bar(dchTimeMin)
title('Decoherence Period Length (min)');
goodFigPrefs

%% Plot start Delay of decoherence period
dchStartDelay = cell2mat(dch.decoherenceStartDelay)./60;
figure();
bar(dchStartDelay)
title('Start Delay of decoherence period(min)');
goodFigPrefs