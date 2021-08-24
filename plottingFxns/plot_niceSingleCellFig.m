function plot_niceSingleCellFig(allCells, i)

addpath(genpath('./plottingFxns'))
 %%
if isempty(i)
    i = randi(size(allCells.spatialFR4,1));
end
% i = 1148 %  G4, WT
% i = 300 %  G2, WT
% i = 3213%  HCNd1, WT
% i = 3484 %  HCNd2, KO
% i = 545 % G3, WT
% i = 285 % G2, WT
% i = 297 % G2, WT
% i = 302 % G2, WT
% i = 481 % G3, WT
% i = 258 %G2, WT
% i = 297 % G2,WT
% i = 689
% i = 566; %G3, WT
%%

singleCellFR4cm = squeeze(allCells.spatialFR4(i,:,:));
singleCellsmoothFR = squeeze(allCells.spatialFRsmooth(i,:,:));
singleCellFR2cm = squeeze(allCells.spatialFR2(i,:,:));
name = allCells.metadata{i,2};
genotype = allCells.metadata{i,4};
% 
ogSpatialBinSize = 2;
spatialBinSize = 4;
numCol2AvgOver = spatialBinSize/ogSpatialBinSize;
singleCellFR = reshape(nanmean(reshape(singleCellFR2cm.',numCol2AvgOver,[])),size(singleCellFR2cm,2)/numCol2AvgOver,[]).';

lickx = allCells.lickX(i).lickx;
lickt = allCells.lickT(i).lickt;
posx = allCells.posX(i).posx;
post = allCells.posT(i).post;
spike_idx = allCells.spike_idx(i);
spike_idx = spike_idx{1};
% FRtime = allCells.FRtime(i).FRtime';
[~,~,lick_idx] = histcounts(lickt,post);

speed = allCells.speed(i).speed;
trial = allCells.trial(i).trial;

dch = allCells.dch(i).dch;

meanSpeed = nan(max(trial),1);
lickAccuracyByTrial = zeros(1,max(trial));
for j = 1:max(trial)
    meanSpeed(j) = nanmean(speed(trial == j));

    trialLicks = lickx(trial(lick_idx) == j);
    goodLicks = sum(trialLicks<25) + sum(trialLicks>max(posx)-25); 
    if trialLicks ~= 0
        lickAccuracyByTrial(j) = goodLicks/numel(trialLicks);
    else
        lickAccuracyByTrial(j) = 0.0;
    end
end

linearFractionalOccupancy = calculate_1D_LFO(posx,post,spatialBinSize,speed);
[I_sec, I_spike] = calculate_1DspatialInformation(singleCellFR4cm,linearFractionalOccupancy);


trialBlockSpatialInformation = nan(max(trial),2);
trialAvgSpeed = nan(max(trial),1);
for j = 1:max(trial)
    trial_posx = posx(trial==j);
    trial_post = post(trial==j);
    trial_speed = speed(trial==j);
    linearFractionalOccupancyBlock = calculate_1D_LFO(trial_posx,trial_post,spatialBinSize,trial_speed);

    [I_sec_block, I_spike_block] = calculate_1DspatialInformation(singleCellFR4cm(j,:),linearFractionalOccupancyBlock);
    trialBlockSpatialInformation(j,:) = [I_sec_block, I_spike_block];
    trialAvgSpeed(j) = nanmean(trial_speed);
end

% Stability Score Curve
stabilityScore = calculateStabilityScore(singleCellFR4cm);

% Correlation Matrix 
%P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
[corrMatrix,pval] = corr(singleCellFR4cm'); 


% Off Diagonal Stability
offdiag = nan(300,1);
for x = 6:295
   offdiag(x) = mean(corrMatrix(x,x-5:x+5)); 
end

% Peakiness
peakiness = peak2rms(singleCellFR4cm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
plotWidth = 1500;
plotHeight = 600;
titleFontSize = 15;
figure('Position',[200 200 plotWidth plotHeight]);
clf;

row = 6;
col = 4;

columnSubplot = [1,col+1,2*col+1,3*col+1,4*col+1,5*col+1];
subplot(row,col,columnSubplot); %raster
scatter(posx(spike_idx),trial(spike_idx),'k.');
set(gca, 'YDir','reverse')
ylim([0 max(trial)+1]);
xlim([0 400]);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',titleFontSize);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Cell %d: %s,%s',i,name,genotype))
xlabel('VR cm')
ylabel('Trial Number')

subplot(row,col,columnSubplot+1); %ratemap
imagesc(singleCellsmoothFR);
colormap('hot')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
% axis off;
set(gca,'Yticklabel',[])
set(gca,'Xticklabel',[])
set(gca,'FontSize',titleFontSize);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title(sprintf('Max FR: %0.1f Hz',max(singleCellsmoothFR,[],'all')));
title('Spatial FR')
colorbar

subplot(row,col,columnSubplot+2); %speed
imagesc(trialAvgSpeed)
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'Yticklabel',[])
set(gca,'Xticklabel',[])
set(gca,'layer','bottom');
set(gca,'FontName','Helvetica');
title('Speed')
box off;
set(gca,'FontSize',titleFontSize);
ylabel('')
colorbar

subplot(row,col,columnSubplot+3) %licks
scatter(posx(lick_idx),trial(lick_idx),'k');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.125]);
set(gca,'Yticklabel',[])
set(gca,'Xticklabel',[])
box off;
set(gca,'layer','bottom');
set(gca,'FontName','Helvetica');
title('Licks')
set(gca,'FontSize',titleFontSize);
ylabel('')
% xlabel('VR cm')
% set(gca,'FontSize',10);

% %%
% subplot(row,col,columnSubplot+2);
% imagesc(dch.idxCellArray);
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis off;
% set(gca,'FontSize',titleFontSize);
% set(gca,'FontName','Helvetica');
% % set(gcf,'Position',[100 100 1000 1000])
% title('UMAP clusters')
% ylabel('Trial Number')
% % [unq, ~, iunq] = unique(dch.idxCellArray);
% % ncol = max(dch.idxCellArray);
% % cb = colorbar;
% % set(gca, 'clim', [0.5 ncol+0.5]);
% % set(cb, 'ticks', 1:1:ncol, 'ticklabels', cellstr(num2str(unq)));
% % set(gca,'TickDir','out');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spatial Firing Rate Map at 4 timepoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
rowStart = 1;
rowEnd = 2;
ylim_range = [0 max(max(mean(smoothdata(singleCellFR2cm(:,:)),2)))];
smoothFactor = 10;

row = 6;
col = 2;
subplot(row,col,rowStart:rowEnd);
plot_lineWithSEM(smoothdata(singleCellFR2cm(1:50,:),2,'movmean',smoothFactor),[])
title('Trial 1-50');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',titleFontSize);
set(gca,'xtick',[]);
ylim(ylim_range);

subplot(row,col,rowStart+col:rowEnd+col)
plot_lineWithSEM(smoothdata(singleCellFR2cm(51:100,:),2,'movmean',smoothFactor),[])
title('Trial 51-100');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',titleFontSize);
set(gca,'xtick',[]);
ylim(ylim_range);

rowPlotNum = 2;
subplot(row,col,rowStart+col*rowPlotNum:rowEnd+col*rowPlotNum)
plot_lineWithSEM(smoothdata(singleCellFR2cm(101:150,:),2,'movmean',smoothFactor),[])
title('Trial 101-150');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',titleFontSize);
set(gca,'xtick',[]);
ylim(ylim_range);

rowPlotNum = 3;
subplot(row,col,rowStart+col*rowPlotNum:rowEnd+col*rowPlotNum)
plot_lineWithSEM(smoothdata(singleCellFR2cm(151:200,:),2,'movmean',smoothFactor),[])
title('Trial 151-200');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',titleFontSize);
set(gca,'xtick',[]);
ylim(ylim_range);

rowPlotNum = 4;
subplot(row,col,rowStart+col*rowPlotNum:rowEnd+col*rowPlotNum)
plot_lineWithSEM(smoothdata(singleCellFR2cm(201:250,:),2,'movmean',smoothFactor),[])
title('Trial 201-250');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',titleFontSize);
set(gca,'xtick',[]);
ylim(ylim_range);

rowPlotNum = 5;
subplot(row,col,rowStart+col*rowPlotNum:rowEnd+col*rowPlotNum)
plot_lineWithSEM(smoothdata(singleCellFR2cm(251:300,:),2,'movmean',smoothFactor),[])
title('Trial 251-300');
ylabel('Hz')
xlabel('VR cm')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',titleFontSize);
x_tick_label = get(gca,'xticklabels');
new_x_tick_label = cellfun(@(x) str2num(x)*2,x_tick_label);
xticklabels(new_x_tick_label)
ylim(ylim_range);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firing Rate, Crest Factor, Spatial Info, Stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smoothFactor = 5;
linewidth = 1;
plotWidth = 600;
plotHeight = 600;

figure('Position',[200 200 plotWidth plotHeight]);
t = tiledlayout(4,1);

nexttile(t)
plot(smooth(mean(singleCellFR2cm,2),smoothFactor),'k','LineWidth',linewidth);
% plot(smooth(mean(singleCellFR2cm,2),'sgolay'),'r');
title('Avg Trial Firing Rate');
ylabel('Hz');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'Xticklabel',[])
set(gca,'FontSize',12);

% subplot(row,col,rowStart+col:rowEnd+col)
nexttile(t)
plot(smooth(peakiness,smoothFactor),'k','LineWidth',linewidth);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
title('Crest Factor')
ylabel('Peak2RMS Ratio')
set(gca,'Xticklabel',[])
set(gca,'FontSize',12);

% rowPlotNum = 2;
% subplot(row,col,rowStart+col*rowPlotNum:rowEnd+col*rowPlotNum)
nexttile(t)
% plot(trialBlockSpatialInformation(:,1),'k')
plot(smooth(trialBlockSpatialInformation(:,1),smoothFactor),'k','LineWidth',linewidth)
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
title('Spatial Information')
ylabel('SI Score')
set(gca,'Xticklabel',[])
set(gca,'FontSize',12);

% rowPlotNum = 3;
% subplot(row,col,rowStart+col*rowPlotNum:rowEnd+col*rowPlotNum)
% nexttile(t)
% % plot(stabilityScore,'k');
% plot(smooth(stabilityScore,smoothFactor),'k','LineWidth',linewidth);
% % plot(smooth(stabilityScore,10,'sgolay'));
% % plot(stabilityScoreCurve);
% % plot(smooth(stabilityScoreCurve,10,'sgolay'));
% % plot(offdiag,'m')
% % set(gca,'TickDir','out');
% % set(gca,'ticklength',[0.015 0.025]);
% % set(gca,'layer','bottom');
% % box off;
% % ylim([0 400]);
% % set(gca,'FontName','Helvetica');
% title('Stability Score')
% ylabel('Stability Score')
% % set(gca,'Xticklabel',[])
% set(gca,'FontSize',12);
% xlabel('Trial Number')


nexttile(t)
% subplot(row,col,rowStart+col*rowPlotNum:rowEnd+col*rowPlotNum)
plot(smooth(lickAccuracyByTrial,smoothFactor),'k','LineWidth',linewidth);
% scatter(trial(lick_idx),posx(lick_idx),'k');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
title('Lick Accuracy')
xlabel('Trial Number')
ylabel('Lick Accuracy')
set(gca,'FontSize',10);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure(3); clf; 
% subplot(row,col,1:18); 
imagesc(fillmissing(corrMatrix,'linear')); colorbar;

% set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Trial by Trial Correlation Matrix'))
ylabel('Trials')
xlabel('Trials')
goodFigPrefs;
colormap('bone');
colorbar;
%
% subplot(row,col,19:24);
% plot(smooth(offdiag,'sgolay'))
% title('Off Diagonal Correlation');

end


