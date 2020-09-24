function plotDecoherencePlots(cells)


pool_decoherence_sessions = poolDecoherenceSessions(cells);
seshNum = numel({pool_decoherence_sessions.seshName});

%% Plot Cluster Groups by Trials Array
% plot_multiDimensionalCellArray({pool_decoherence_sessions.idxCellArray}')
% goodFigPrefs
% colormap jet
% title('Cluster Groups by Trials Array');

%% Highlight the Decoherence Period for Trials

highlightedDecoherenceIndx = cell(seshNum,1);
for i = 1:seshNum
    iCA = pool_decoherence_sessions(i).idxCellArray;
    iDI = pool_decoherence_sessions(i).decoherenceIdx;
    highlightedDecoherenceIndx{i} = iCA;
    highlightedDecoherenceIndx{i}(iDI) = 10;
end
plot_multiDimensionalCellArray(highlightedDecoherenceIndx)
goodFigPrefs
title('Decoherence Period');
colormap jet  

%%  Plot Cluster Groups by Time Array
plot_multiDimensionalCellArray({pool_decoherence_sessions.idxClusterArray}')
goodFigPrefs

title('Cluster Groups by Time Array');
colormap jet
yt = get(gca, 'YTick'); 
set(gca, 'YTick',yt, 'YTickLabel',round(yt*(1/pool_decoherence_sessions(1).Fs)/60))
ylabel('Minutes');

%% Highlight the Decoherence Period for Time
highlightedDecoherenceIndx = cell(seshNum,1);
for i = 1:seshNum
    iCA = pool_decoherence_sessions(i).idxClusterArray;
    iDI = pool_decoherence_sessions(i).decoherenceTimeIdx;
    iDI = iDI * pool_decoherence_sessions(i).Fs;
    highlightedDecoherenceIndx{i} = iCA;
    highlightedDecoherenceIndx{i}(min(iDI):max(iDI)) = 15;
end

plot_multiDimensionalCellArray(highlightedDecoherenceIndx)
goodFigPrefs
title('Decoherence Period');
colormap jet
yt = get(gca, 'YTick'); 
set(gca, 'YTick',yt, 'YTickLabel',round(yt*(1/pool_decoherence_sessions(1).Fs)/60))
ylabel('Minutes');
xlabel('Session');
%% Plot how many groups UMAP identified per session
figure()
histogram(cellfun(@max,{pool_decoherence_sessions.idxCellArray}'),'FaceColor','black')
goodFigPrefs
title('Number of UMAP identified groups');

%% Plot number of trials
figure();
histogram(cellfun(@numel,{pool_decoherence_sessions.decoherenceIdx}'),10,'FaceColor','black')
goodFigPrefs
title('Decoherence Period Length (trials)');

%% Plot length of decoherence period
dchTimeMin = cell2mat({pool_decoherence_sessions.decoherenceTime}')./60;
figure();
histogram(dchTimeMin,10,'FaceColor','black')
title('Decoherence Period Length (min)');
goodFigPrefs
xlabel('Min')
ylabel('Number of Sessions')

%% Plot start Delay of decoherence period
dchStartDelay = cell2mat({pool_decoherence_sessions.decoherenceStartDelay}')./60;
figure();
histogram(dchStartDelay,10,'FaceColor','black')
title('Start Delay of decoherence period(min)');
goodFigPrefs
xlabel('Min')
ylabel('Number of Sessions')


%% Plot Avg FR before and after Dch Period
dchStartDelay = cell2mat({pool_decoherence_sessions.decoherenceStartDelay}')./60;
figure();
histogram(dchStartDelay,10,'FaceColor','black')
title('Start Delay of decoherence period(min)');
goodFigPrefs
xlabel('Min')
ylabel('Number of Sessions')


end
