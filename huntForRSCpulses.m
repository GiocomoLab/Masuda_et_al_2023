%%load data
dataset = load('/Users/KeiMasuda/Desktop/fkm_analysis/AA_190830_1_0927_baseline+ketamine.mat');
% dataset = load('/Users/KeiMasuda/Desktop/fkm_analysis/AA_190709_2_0812_baseline+ketamine.mat');
% dataset = load('/Users/KeiMasuda/Desktop/fkm_analysis/AA_190709_4_RSC_baseline+ketamine.mat');
%% compute metrics (Depth)
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(dataset.sp.temps, dataset.sp.winv, dataset.sp.ycoords, dataset.sp.spikeTemplates, dataset.sp.tempScalingAmps);
depth=zeros(length(dataset.sp.cgs),1);
for iC=1:length(dataset.sp.cgs)
    depth(iC)=mean(spikeDepths(dataset.sp.clu==dataset.sp.cids(iC)));
end
dataset.sp.spikeAmps = spikeAmps;
dataset.spikeDepths = spikeDepths;
%% generate a spike matrix that is n_depth bins by n_timebins
ketamineInjxTime = dataset.post(find(dataset.trial==30,1));
valid_idx = ismember(dataset.sp.clu,dataset.sp.cids(dataset.sp.cgs==2));
% 
preKetamineIndx = dataset.sp.st>ketamineInjxTime;
postKetamineIndx = dataset.sp.st<ketamineInjxTime;

%%
valid_idx = ismember(dataset.sp.clu,dataset.sp.cids(dataset.sp.cgs==2));
valid_idx(preKetamineIndx) = 0;
bins = 0:40:3840;

discrete_depth = discretize(spikeDepths(valid_idx),bins); %depth bin for each spike
spike_times = dataset.sp.st(valid_idx);
spike_times = double(spike_times);

time_bins = 0.002;
discrete_time = round(spike_times/time_bins)+1;
spikeMat = zeros(numel(bins),(ceil(max(spike_times))/time_bins)+1,'single');
% populate matrix
for iS=1:numel(spike_times)
    time_idx=discrete_time(iS);
    d_idx = discrete_depth(iS);
    
    if ~isnan(d_idx)
        spikeMat(d_idx,time_idx)=spikeMat(d_idx,time_idx)+1;
    end
    
end
spikeMat_preKet = spikeMat;

% calculate power spectrum for each depth
[PxxSpikes,FSpikes] = pwelch(spikeMat_preKet',[],[],[],1/time_bins);
% get power in delta range, normalized by other frequencies
% delta_range=[4 12];
delta_range=[1 3];
delta_idx = FSpikes>delta_range(1) & FSpikes<=delta_range(2);
rest_idx = ~delta_idx;
deltaPower = mean(PxxSpikes(delta_idx,:));
restPower = mean(PxxSpikes(rest_idx,:));
deltaPowerN = deltaPower./restPower;
%find depth with strongest delta power
[a,maxChan_spikes]=max(deltaPowerN);
%bandpass that vector around delta band
bp_spikes = bandpass(spikeMat(maxChan_spikes,:),delta_range,1/time_bins);
%find peaks
[pks,locs_spikes] = findpeaks(bp_spikes);

% extract snippets around peaks in all channels
snpsSpikes=extract_snps(spikeMat,locs_spikes,'win',[-100 100]);
tvec_spikes = [-100:100]*time_bins;
aa_spikes=squeeze(mean(snpsSpikes,3));

% some basic plotting
aa_spikesNorm=bsxfun(@rdivide,aa_spikes,sum(aa_spikes,2));

figure(1)
subplot(1,2,1)
imagesc(flipud(aa_spikesNorm))
set(gca,'YTick',linspace(1,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,numel(tvec_spikes),5),'XTickLabel',linspace(min(tvec_spikes),max(tvec_spikes),5))
%yline(385-highest_channel,'k');
vline(numel(tvec_spikes)/2+.5)
title('pre ketamine triggered spikes')
[~,max_loc]=max(aa_spikes,[],2);
hold on
plot((max_loc(end:-1:1)),1:numel(max_loc),'ro')

set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;

set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])

% ===================================================
% calculate power spectrum for each depth for POST KET
% ===================================================

% Calculate for postKetamineIndx
valid_idx = ismember(dataset.sp.clu,dataset.sp.cids(dataset.sp.cgs==2));
valid_idx(postKetamineIndx) = 0;
bins = 0:40:3840;

discrete_depth = discretize(spikeDepths(valid_idx),bins); %depth bin for each spike
spike_times = dataset.sp.st(valid_idx);
spike_times = spike_times - ketamineInjxTime;
spike_times = double(spike_times);

time_bins = 0.002;
discrete_time = round(spike_times/time_bins)+1;
spikeMat = zeros(numel(bins),(ceil(max(spike_times))/time_bins)+1,'single');
% populate matrix
for iS=1:numel(spike_times)

    time_idx=discrete_time(iS);
    d_idx = discrete_depth(iS);
    
    if ~isnan(d_idx)
        spikeMat(d_idx,time_idx)=spikeMat(d_idx,time_idx)+1;
    end
    
end
spikeMat_postKet = spikeMat;



[PxxSpikes,FSpikes] = pwelch(spikeMat_postKet',[],[],[],1/time_bins);
% get power in delta range, normalized by other frequencies
% delta_range=[4 12];
delta_range=[1 3];
delta_idx = FSpikes>delta_range(1) & FSpikes<=delta_range(2);
rest_idx = ~delta_idx;
deltaPower = mean(PxxSpikes(delta_idx,:));
restPower = mean(PxxSpikes(rest_idx,:));
deltaPowerN = deltaPower./restPower;
%find depth with strongest delta power
[a,maxChan_spikes]=max(deltaPowerN);
%bandpass that vector around delta band
bp_spikes = bandpass(spikeMat(maxChan_spikes,:),delta_range,1/time_bins);
%find peaks
[pks,locs_spikes] = findpeaks(bp_spikes);

% extract snippets around peaks in all channels
snpsSpikes=extract_snps(spikeMat,locs_spikes,'win',[-100 100]);
tvec_spikes = [-100:100]*time_bins;
aa_spikes=squeeze(mean(snpsSpikes,3));

% some basic plotting
aa_spikesNorm=bsxfun(@rdivide,aa_spikes,sum(aa_spikes,2));

subplot(1,2,2)
imagesc(flipud(aa_spikesNorm))
set(gca,'YTick',linspace(1,size(spikeMat,1),10),'YTickLabel',round(linspace(max(bins),min(bins),10)))
set(gca,'XTick',linspace(1,numel(tvec_spikes),5),'XTickLabel',linspace(min(tvec_spikes),max(tvec_spikes),5))
%yline(385-highest_channel,'k');
vline(numel(tvec_spikes)/2+.5)
title('post ketamine triggered spikes')
[~,max_loc]=max(aa_spikes,[],2);
hold on
plot((max_loc(end:-1:1)),1:numel(max_loc),'ro')


set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;

set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])


%%
Fs = 1/time_bins;
depthRange = 60:80;
% depthRange = 1:size(spikeMat,1);


timeRangeCircaKetamineMin= 7;
timeRangeCircaKetamineIndx = timeRangeCircaKetamineMin*60*Fs;

close(figure(2)); figure(2);
timeSize = size(spikeMat_preKet,2);
[pxx,w] = pwelch(nanmean(spikeMat_preKet(depthRange,:),1)',[],[],[],Fs,'onesided','psd');
smoothing = size(pxx,1)/100;
% smoothing = 100;

x = w/pi;
y = smooth(10*log10(pxx),smoothing,'sgolay');
plot(x(x<10),y(x<10),'lineWidth',2,'displayName','Pre Ketamine');

hold on;

pwelchData = nanmean(spikeMat_postKet(depthRange,:),1);
[pxx,w] = pwelch(pwelchData,[],[],[],Fs,'onesided','psd');
x = w/pi;
y = smooth(10*log10(pxx),smoothing,'sgolay');
plot(x(x<10),y(x<10),'lineWidth',2,'displayName','Post Ketamine');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',15);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title('Periodogram Pre vs Post Ketamine')
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
legend
% ylim([-50 -25])