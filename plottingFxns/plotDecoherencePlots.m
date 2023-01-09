 function plotDecoherencePlots(cells)

addpath(genpath(pwd))
pool_decoherence_sessions = poolDecoherenceSessions(cells);
seshNum = numel({pool_decoherence_sessions.seshName});

%% Plot Cluster Groups by Trials Array
% plot_multiDimensionalCellArray({pool_decoherence_sessions.idxCellArray}')
% goodFigPrefs
% colormap jet
% title('Cluster Groups by Trials Array');

%% Highlight the Decoherence Period for Trials
figure()
highlightedDecoherenceIndx = cell(seshNum,1);
for i = 1:seshNum
    try
        iCA = pool_decoherence_sessions(i).idxCellArray;
        iDI = pool_decoherence_sessions(i).decoherenceIdx;
        highlightedDecoherenceIndx{i} = iCA;
        highlightedDecoherenceIndx{i}(iDI) = 20;
    catch
        fprintf('\nFailed')
    end
end
plot_multiDimensionalCellArray(highlightedDecoherenceIndx)
goodFigPrefs
title('Decoherence Period');
colormap jet  
ylabel('Trials');
xlabel('Session');
%%  Plot Cluster Groups by Time Array
% plot_multiDimensionalCellArray({pool_decoherence_sessions.idxClusterArray}')
% goodFigPrefs
% 
% title('Cluster Groups by Time Array');
% colormap jet
% yt = get(gca, 'YTick'); 
% set(gca, 'YTick',yt, 'YTickLabel',round(yt*(1/pool_decoherence_sessions(1).Fs)/60))
% ylabel('Minutes');

%% Highlight the Decoherence Period for Time
figure()
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
try
figure()
histogram(cellfun(@max,{pool_decoherence_sessions.idxCellArray}'),'FaceColor','black')
goodFigPrefs
title('Number of UMAP identified groups');
ylabel('Number of Sessions');
xlabel('Num of UMAP Clusters');
end
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
histogram(dchStartDelay,15,'FaceColor','black')
title('Start Delay of decoherence period(min)');
goodFigPrefs
xlabel('Min')
ylabel('Number of Sessions')

%% Plot Stats and Bar Graphs comparing before and after dch period
plot_STATS_before_and_after_dchPeriod(cells)

%% a plot that looks for spatial autocorrelation in the deocherence period
plot_autocorrelation_before_and_after_dchPeriod(cells)

end
