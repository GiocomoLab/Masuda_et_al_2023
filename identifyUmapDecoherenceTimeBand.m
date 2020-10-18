function [dch_idx, time_idx, dchTimeSec,dchStartDelaySec,Fs, trialClust,trialBasedClustIdentifiers]...
    = identifyUmapDecoherenceTimeBand(clusterIdentifiers, cells, ds_factor)
% Threshold Variables
% return cluster of dchTime is longer than 1 min and 
% if it starts within 10 min of the ketamine injection
dchTimeThreshold_min = 2;% min length of decoherence period
dchLengthMaxTime_min = 60; %max length of decoherence period 
dchStartDelay_min = 15 ; % maximum start delay of decoherence period

%% IDs the decoherence band from the results of UMAP being run on a
% time-binned FR across all cells
trials = downsample(cells.trial(1).trial,ds_factor);
post = downsample(cells.posT(1).post,ds_factor);
Fs = 1/(post(2)-post(1)); %sampling rate is 50hz

ketamineTrial = 100; 
%%
% Create a vector of postKet trial length and assigns a cluster value to each trial based
% on the mode cluster value for the trial duration
minTrial = min(trials);
maxTrial = max(trials);
trialClust = nan(maxTrial-minTrial,1);
for k = minTrial:maxTrial
    trialClust(k) = mode(clusterIdentifiers(trials==k));
end

trialBasedClustIdentifiers = nan(size(trials));
for idx = 1:numel(trials)
    trialBasedClustIdentifiers(idx) = trialClust(trials(idx));
end
%%
% Find time indices of the different smoothed clusters
clusterBoundaries = diff(trialClust') ~= 0;
clstBegin = find([true,clusterBoundaries])';  % beginning indices
clstEnd = find([clusterBoundaries,true])';    % ending indices
clstLen = 1 + clstEnd - clstBegin;            % cluster length
clst = [clstBegin, clstEnd, clstLen];

% return the indices of the cluster of interest
dch_idx = []; time_idx=[]; dchTimeSec=[]; dchStartDelaySec=[]; 

for j = 1:size(clst,1)
   row = clst(j,:);
   dch_trialIdx = row(1):row(2);

   if dch_trialIdx(1) > ketamineTrial

       % Find Decoherence period length
       dchTimeIndxStart = find(trials==row(1),1,'first');
       dchTimeIndxEnd = find(trials==row(2),1,'last');
       dchTimeIndx = dchTimeIndxStart:dchTimeIndxEnd;
       dchTimeSec_temp = numel(dchTimeIndx)/Fs;    

       % Find Start Delay - this is the time the cluster started;
       % this works when the clusters being passed into this already had
       % the first 100 pre-ketamine trials removed
       ketamineTimeIndx = find(trials>ketamineTrial,1,'first');
       dchStartDelaySec_temp = (dchTimeIndxStart-ketamineTimeIndx)/Fs;


       % return first cluster of dchTime that is longer than 5 min and 
       % if it starts within 15 min of the ketamine injection
       if dchTimeSec_temp>dchTimeThreshold_min*60 && dchStartDelaySec_temp<dchStartDelay_min*60
           % make sure decoherence period is less than 1 hour
           if dchTimeSec_temp<dchLengthMaxTime_min*60
               fprintf('Found dch period of %.2f min\n',dchTimeSec_temp/60);
               dch_idx = dch_trialIdx;
               %convert dechorence index values into seconds
               time_idx = dchTimeIndx./Fs;
               dchTimeSec=dchTimeSec_temp;
               dchStartDelaySec=dchStartDelaySec_temp;
               fprintf('Found Decoherence Band: Cluster %i, Length: %3.1fmin\n',j,dchTimeSec/60)
               break
           end
       end
   end
end     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot 1 example raster
close all;
dchRange = dch_idx;
plotWidth = 160;
plotHeight = 500;

h = figure('Position',[100 100 plotWidth plotHeight]); hold on;
i = randi(size(cells.spatialFR2,1),1);

name = cells.metadata{i,2};
genotype = cells.metadata{i,4};
sessionDate = cells.metadata{i,3};

posx = cells.posX(i).posx;
spike_idx = cells.spike_idx(i);
spike_idx = spike_idx{1};
% FRtime = allCells.FRtime(i).FRtime';

trial = cells.trial(i).trial;


clf;
scatter(posx(spike_idx),trial(spike_idx),0.5,'k.'); hold on;

if ~isempty(dchRange)
    x = posx(spike_idx);
    y = trial(spike_idx);
    x = x(y>min(dchRange) & y<max(dchRange));
    y = y(y>min(dchRange) & y<max(dchRange));
    scatter(x,y,0.5,'r.');
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

%%
figure('Position',[500 500 plotWidth plotHeight]); clf;
imagesc(trialClust);

end

