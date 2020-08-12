function [dch_idx, time_idx, dchTimeSec,dchStartDelaySec, Fs, trialClust] = identifyUmapDecoherenceTimeBand(clusterIdentifiers, cells, ds_factor)
% IDs the decoherence band from the results of UMAP being run on a
% time-binned FR across all cells

trials = downsample(cells.trial(1).trial,ds_factor);
post = downsample(cells.posT(1).post,ds_factor);
Fs = 1/(post(2)-post(1)); %sampling rate is 50hz


% smooth clusterIdentifiers with a moving median(window of 1sec/10samples
[smth_clusterIdentifiers,window] = smoothdata(clusterIdentifiers,'movmedian',50);
minTrial = min(trials);
maxTrial = max(trials);
trialClust = nan(maxTrial-minTrial,1);
for k = minTrial:maxTrial
    trialClust(k) = median(smth_clusterIdentifiers(trials==k));
end

figure(1)
scatter(trials,smth_clusterIdentifiers,100,smth_clusterIdentifiers)
colormap('jet');


% Find index trial cluster regions
clusterBoundaries = diff(smth_clusterIdentifiers) ~= 0;
clstBegin = find([true,clusterBoundaries])';  % beginning indices
clstEnd = find([clusterBoundaries,true])';    % ending indices
clstLen = 1 + clstEnd - clstBegin;            % cluster length

clst = [clstBegin, clstEnd, clstLen];

% Pull out first cluster after trial 100
postKetTrialsIdx = find(trials > 100,1,'first');
% Pull out first cluster after trial 100 
postKetClst = clst(clst(:,1) >= postKetTrialsIdx,:);


% return the indices of the cluster of interest
dch_idx = []; time_idx=[]; dchTimeSec=[]; dchStartDelaySec=[]; 
try
   for j = 1:size(postKetClst,1)
       row = postKetClst(j,:);
       dch_idx_temp = row(1):row(2);
       
       % Find Decoherence period length
       dchTimeSec_temp = numel(dch_idx_temp)/Fs;
%        dchTimeMin = dchTimeSec_temp/60
        
       % Find Start Delay
       idxTrialStart = row(1);
       idxKetamineStart = postKetTrialsIdx;
       dchStartDelaySec_temp = (idxTrialStart-idxKetamineStart)/Fs;
%        dchStartDelayMin = dchStartDelaySec_temp/60;
       
       % return cluster of dchTime is longer than 5 min and 
       % if it starts within 2 min of the ketamine injection
       if dchTimeSec_temp>5*60 && dchStartDelaySec_temp<5*60
           fprintf('Found dch period of %.2f min\n',dchTimeSec_temp/60);
           dch_idx = min(trials(dch_idx_temp)):max(trials(dch_idx_temp));
           time_idx = dch_idx_temp./Fs;
           dchTimeSec=dchTimeSec_temp;
           dchStartDelaySec=dchStartDelaySec_temp;
           break
       end
   end
   
catch
   dch_idx = []; dchTimeSec=[]; dchStartDelaySec=[];
end


end

% Keep only cluster of length 10+trials
% clstIdx = find(clstLen > 3);
% clst = [clstBegin(clstIdx), clstEnd(clstIdx), clstLen(clstIdx)];
% uniqueVal = unique(clst(:,4))

