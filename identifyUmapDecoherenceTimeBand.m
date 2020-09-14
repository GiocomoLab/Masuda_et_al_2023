function [dch_idx, time_idx, dchTimeSec,dchStartDelaySec, smth_clusterIdentifiers_all,Fs, trialClust]...
    = identifyUmapDecoherenceTimeBand(clusterIdentifiers, cells, ds_factor)
% Threshold Variables
% return cluster of dchTime is longer than 1 min and 
% if it starts within 10 min of the ketamine injection
dchTimeThreshold_min = 5;% min length of decoherence period
dchLengthThreshold_min = 60; %max length of decoherence period 
dchStartDelay_min = 15 ; % maximum start delay of decoherence period

%% IDs the decoherence band from the results of UMAP being run on a
% time-binned FR across all cells
trials = downsample(cells.trial(1).trial,ds_factor);
post = downsample(cells.posT(1).post,ds_factor);
Fs = 1/(post(2)-post(1)); %sampling rate is 50hz

ketamineTrial = 100; 
% Pull out first idx after trial 100 (ketamine admin after trial 100)
postKetTrialsIdx = find(trials > ketamineTrial,1,'first');
% Pull out clusters after ketamine injection (trial 100)
postKetClusterIdentifiers = clusterIdentifiers(1,postKetTrialsIdx:end);
% Pull out trial indices after ketamine injection
postKetTrials = trials(postKetTrialsIdx:end);

%% smooth clusterIdentifiers with a moving median(window is the downsampleing factor)
clusterSmoothingConstant = ds_factor;
smth_clusterIdentifiers = round(smoothdata(postKetClusterIdentifiers,'movmedian',clusterSmoothingConstant));
smth_clusterIdentifiers_all = round(smoothdata(clusterIdentifiers,'movmedian',clusterSmoothingConstant));
%
% figure()
% scatter(postKetTrials,smth_clusterIdentifiers,100,postKetClusterIdentifiers)
% colormap('jet');
% goodFigPrefs();
% title('Smoothed IDed Clusters after trial 100');
% xlabel('Trial');
% ylabel('UMAP Cluster Group');
% colorbar();
%%

% Create a vector of postKet trial length and assigns a cluster value to each trial based
% on the median cluster value for the trial duration
minTrial = min(postKetTrials);
maxTrial = max(trials);
trialClust = nan(maxTrial-minTrial,1);
for k = minTrial:maxTrial
    trialClust(k) = median(smth_clusterIdentifiers(postKetTrials==k));
end

%%
% Find time indices of the different smoothed clusters
clusterBoundaries = diff(smth_clusterIdentifiers) ~= 0;
clstBegin = find([true,clusterBoundaries])';  % beginning indices
clstEnd = find([clusterBoundaries,true])';    % ending indices
clstLen = 1 + clstEnd - clstBegin;            % cluster length
clst = [clstBegin, clstEnd, clstLen];

% return the indices of the cluster of interest
dch_idx = []; time_idx=[]; dchTimeSec=[]; dchStartDelaySec=[]; 
try
   for j = 1:size(clst,1)
       row = clst(j,:);
       dch_idx_temp = row(1):row(2);
       
       % Find Decoherence period length
       dchTimeSec_temp = numel(dch_idx_temp)/Fs;
%        dchTimeMin = dchTimeSec_temp/60
        
       % Find Start Delay
       idxTrialStart = row(1);
       idxKetamineStart = postKetTrialsIdx;
       dchStartDelaySec_temp = (idxTrialStart-idxKetamineStart)/Fs;
%        dchStartDelayMin = dchStartDelaySec_temp/60;
       
       % return first cluster of dchTime that is longer than 5 min and 
       % if it starts within 15 min of the ketamine injection
       if dchTimeSec_temp>dchTimeThreshold_min*60 && dchStartDelaySec_temp<dchStartDelay_min*60
           % make sure cluster is less than 1 hour
           if dchTimeSec_temp<dchLengthThreshold_min*60
               fprintf('Found dch period of %.2f min\n',dchTimeSec_temp/60);
               dch_idx = min(postKetTrials(dch_idx_temp)):max(postKetTrials(dch_idx_temp));
               time_idx = dch_idx_temp./Fs;
               dchTimeSec=dchTimeSec_temp;
               dchStartDelaySec=dchStartDelaySec_temp;
               break
           end
       end
       fprintf('Band %i does not qualify as a decoherence band\n',j)
   end
   
catch
   dch_idx = []; dchTimeSec=[]; dchStartDelaySec=[];
end

%% plot 1 example raster
close all;
dchRange = dch_idx;
plotWidth = 160;
plotHeight = 500;

h = figure('Position',[100 100 plotWidth plotHeight]); hold on;
i = randi(size(cells.spatialFR10,1),1);


singleCellFR2cm = squeeze(cells.spatialFR2(i,:,:));

name = cells.metadata{i,2};
genotype = cells.metadata{i,4};
sessionDate = cells.metadata{i,3};
% 
ogSpatialBinSize = 2;
spatialBinSize = 10;
numCol2AvgOver = spatialBinSize/ogSpatialBinSize;
singleCellFR = reshape(nanmean(reshape(singleCellFR2cm.',numCol2AvgOver,[])),size(singleCellFR2cm,2)/numCol2AvgOver,[]).';

posx = cells.posX(i).posx;
spike_idx = cells.spike_idx(i);
spike_idx = spike_idx{1};
% FRtime = allCells.FRtime(i).FRtime';

trial = cells.trial(i).trial;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

scatter(posx(spike_idx),trial(spike_idx),1,'k.'); hold on;

if ~isempty(dchRange)
    x = posx(spike_idx);
    y = trial(spike_idx);
    x = x(y>min(dchRange) & y<max(dchRange));
    y = y(y>min(dchRange) & y<max(dchRange));
    scatter(x,y,1,'r.');
end
colormap('default')

set(gca, 'YDir','reverse')
ylim([0 max(trial)+1]);
xlim([0 400]);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.01 0.025]);   
set(gca,'layer','bottom');

set(gca,'FontName','Helvetica');
box off;

set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Cell %d: %s,%s',i,name,genotype))


end

