function plot_niceSingleCellFig(allCells)

addpath(genpath('./plottingFxns'))
%%
i = randi(size(allCells.spatialFR10,1));
% i = 1148 G4, WT
% i = 300  G2, WT
% i = 3213 HCNd1, WT
% i = 3484 HCNd2, KO
% i = 545 G3, WT
% i = 285 G2, WT
% i = 297 G2, WT
% i = 302 G2, WT
% i = 481 G3, WT
%

singleCellFR10cm = squeeze(allCells.spatialFR10(i,:,:));
singleCellFR2cm = squeeze(allCells.spatialFR2(i,:,:));
name = allCells.metadata{i,2};
genotype = allCells.metadata{i,4};
% 
ogSpatialBinSize = 2;
spatialBinSize = 10;
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
[I_sec, I_spike] = calculate_1DspatialInformation(singleCellFR10cm,linearFractionalOccupancy);


trialBlockSpatialInformation = nan(max(trial),2);
for j = 1:max(trial)
    trial_posx = posx(trial==j);
    trial_post = post(trial==j);
    trial_speed = speed(trial==j);
    linearFractionalOccupancyBlock = calculate_1D_LFO(trial_posx,trial_post,spatialBinSize,trial_speed);

    [I_sec_block, I_spike_block] = calculate_1DspatialInformation(singleCellFR10cm(j,:),linearFractionalOccupancyBlock);
    trialBlockSpatialInformation(j,:) = [I_sec_block, I_spike_block];
end

% Stability Score Curve
stabilityScore = calculateStabilityScore(singleCellFR10cm);

% Correlation Matrix 
%P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
[corrMatrix,pval] = corr(singleCellFR10cm'); 


% Off Diagonal Stability
offdiag = nan(300,1);
for x = 6:295
   offdiag(x) = mean(corrMatrix(x,x-5:x+5)); 
end

% Peakiness
peakiness = peak2rms(singleCellFR10cm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
figure(1);
clf;

row = 5;
col = 5;

columnSubplot = [1,6,11,16,21];
subplot(row,col,columnSubplot)
scatter(posx(spike_idx),trial(spike_idx),'k.');
colormap('default')
set(gca, 'YDir','reverse')
ylim([0 max(trial)+1]);
xlim([0 400]);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Cell %d: %s,%s',i,name,genotype))
xlabel('VR cm')
ylabel('Trial Number')

subplot(row,col,columnSubplot+1)
imagesc(singleCellFR);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis off;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title('Firing Rate')
xlabel('VR cm')
ylabel('Trial Number')

subplot(row,col,columnSubplot+2)
imagesc(dch.idxCellArray);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis off;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title('UMAP clusters')
ylabel('Trial Number')
[unq, ~, iunq] = unique(dch.idxCellArray);
ncol = max(dch.idxCellArray);
cb = colorbar;
set(gca, 'clim', [0.5 ncol+0.5]);
set(cb, 'ticks', 1:1:ncol, 'ticklabels', cellstr(num2str(unq)));
set(gca,'TickDir','out');
%%%%%%%
smoothFactor = 10;

rowStart = 4;
rowEnd = 5;
subplot(row,col,rowStart:rowEnd);

plot(smooth(mean(singleCellFR2cm,2),smoothFactor),'k');
% plot(smooth(mean(singleCellFR2cm,2),'sgolay'),'r');
title('Avg Trial Firing Rate');
ylabel('Hz');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'Xticklabel',[])
set(gca,'FontSize',15);

subplot(row,col,rowStart+col:rowEnd+col)
plot(smooth(peakiness,smoothFactor),'k');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
title('Crest Factor')
ylabel('Peak2RMS Ratio')
set(gca,'Xticklabel',[])
set(gca,'FontSize',15);

rowPlotNum = 2;
subplot(row,col,rowStart+col*rowPlotNum:rowEnd+col*rowPlotNum)
plot(smooth(trialBlockSpatialInformation(:,2),smoothFactor),'k')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
title(sprintf('Spatial Information = Isec:%0.3f; Ispike:%0.3f', I_sec, I_spike))
ylabel('SI Score')
set(gca,'Xticklabel',[])
set(gca,'FontSize',15);

rowPlotNum = 3;
subplot(row,col,rowStart+col*rowPlotNum:rowEnd+col*rowPlotNum)
plot(smooth(stabilityScore,smoothFactor),'k');
% plot(smooth(stabilityScore,10,'sgolay'));
% plot(stabilityScoreCurve);
% plot(smooth(stabilityScoreCurve,10,'sgolay'));
% plot(offdiag,'m')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
title('Stability Score')
ylabel('Stability Score')
set(gca,'Xticklabel',[])
set(gca,'FontSize',15);

rowPlotNum = 4;
subplot(row,col,rowStart+col*rowPlotNum:rowEnd+col*rowPlotNum)
% bar(lickAccuracyByTrial);
scatter(trial(lick_idx),posx(lick_idx),'k');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
ylim([0 400]);
set(gca,'FontName','Helvetica');
title('Licks')
xlabel('Trial Number')
ylabel('VR cm')
set(gca,'FontSize',15);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial Firing Rate Map at 4 timepoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf; 
row = 4;
col = 6;
subplot(row,col,1:6);
plot_lineWithSEM(smoothdata(singleCellFR2cm(1:50,:),2,'movmean',smoothFactor),[])
title('Trial 1-50');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',20);
set(gca,'xtick',[]);

subplot(row,col,7:12);
plot_lineWithSEM(smoothdata(singleCellFR2cm(51:100,:),2,'movmean',smoothFactor),[])
title('Trial 51-100');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',20);
set(gca,'xtick',[]);

subplot(row,col,13:18);
plot_lineWithSEM(smoothdata(singleCellFR2cm(101:115,:),2,'movmean',smoothFactor),[])
title('Trial 101-115');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',20);
set(gca,'xtick',[]);
% ylim([0 20])

subplot(row,col,19:24);
plot_lineWithSEM(smoothdata(singleCellFR2cm(116:300,:),2,'movmean',smoothFactor),[])
title('Trial 116-300');
ylabel('Hz')
xlabel('VR cm')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf; 
% subplot(row,col,1:18); 
imagesc(fillmissing(corrMatrix,'linear')); colorbar;
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box on;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica'); 
% set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Trial by Trial Correlation Matrix'))
ylabel('Trials')
xlabel('Trials')

%
% subplot(row,col,19:24);
% plot(smooth(offdiag,'sgolay'))
% title('Off Diagonal Correlation');

end


