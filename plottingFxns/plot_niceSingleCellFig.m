function plot_niceSingleCellFig(allCells)

%%
i = randi(size(allCells.spatialFR10,1));
% i = 1148
% i = 300
% i = 3213
% i = 3484
% i = 545

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
testCells = singleCellFR10cm;
numCells = size(testCells,1);
trialNum = size(testCells,2);
spatialBins = size(testCells,3);
flatFR = reshape(permute(testCells,[3 1 2]), [numCells*spatialBins, trialNum]);
flatFR = normalize(flatFR,1)';

%P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
corrMatrix = corr(fillmissing(flatFR,'linear')); 

% Off Diagonal Stability
offdiag = nan(300,1);
for x = 6:295
   offdiag(x) = mean(corrMatrix(x,x-5:x+5)); 
end

% Peakiness
peakiness = peak2rms(singleCellFR10cm');

%%%%%%%%%%%%%%%%
figure(1);
row = 5;
col = 6;

subplot(row,col,[1,7,13,19,25])
scatter(posx(spike_idx),trial(spike_idx),'k.');
colormap('default')
set(gca, 'YDir','reverse')
ylim([0 max(trial)+1]);
xlim([0 400]);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
% set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Cell %d ? %s,%s',i,name,genotype))
xlabel('VR cm')
ylabel('Trial Number')

subplot(row,col,[2,8,14,20,26])
imagesc(singleCellFR);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
% set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title('Firing Rate')
xlabel('VR cm')
ylabel('Trial Number')



% subplot(row,col,3:6)
% plot(smooth(meanSpeed,'sgolay'));
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% % set(gca,'FontSize',20);
% set(gca,'FontName','Helvetica');
% % set(gcf,'Position',[100 100 1000 1000])
% title('Speed')
% ylabel('cm/s')
%
subplot(row,col,3:6)
plot(smooth(peakiness,'sgolay'));
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
title('Peakiness')
ylabel('Ratio of largest absolute to RMS')


subplot(row,col,9:12); 
plot(smooth(mean(singleCellFR2cm,2),'sgolay'),'r');
title('Avg Trial Firing Rate');
ylabel('Hz');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;

subplot(row,col,15:18);
stabilityScoreCurve = squeeze(allCells.stabilityScoreCurve(i,:));
% plot(stabilityScore);
% plot(smooth(stabilityScore,10,'sgolay'));
plot(smooth(stabilityScoreCurve,10,'sgolay'));
% plot(offdiag,'m')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
% set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title('Stability Score')
ylabel('Stability Score')


subplot(row,col,21:24);
plot(smooth(trialBlockSpatialInformation(:,2),'sgolay'),'k')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
title(sprintf('Spatial Information = Isec:%0.3f; Ispike:%0.3f', I_sec, I_spike))
% xlabel('Trial Number')
ylabel('SI Score')
% legend('I sec', 'I spike');
%


subplot(row,col,27:30);
% bar(lickAccuracyByTrial);
scatter(trial(lick_idx),posx(lick_idx));
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
ylim([0 400]);
set(gca,'FontName','Helvetica');
title('Licks')
xlabel('Trial Number')
ylabel('VR cm')

%%%%%%%%%%%%%%%%%%
figure(2); 
row = 4;
col = 6;
subplot(row,col,1:6);
plot(mean(singleCellFR2cm(1:50,:),1))
title('Trial 1-50');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',15);

subplot(row,col,7:12);
plot(mean(singleCellFR2cm(51:100,:),1))
title('Trial 51-100');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',15);

subplot(row,col,13:18);
plot(mean(singleCellFR2cm(101:125,:),1))
title('Trial 101-125');
ylabel('Hz')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',15);

subplot(row,col,19:24);
plot(mean(singleCellFR2cm(126:300,:),1))
title('Trial 126-300');
ylabel('Hz')
xlabel('VR cm')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',15);

%


figure(3)
subplot(row,col,1:18); 
imagesc(corrMatrix); colorbar;
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box on;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica'); 
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Trial by Trial Correlation Matrix'))
%
subplot(row,col,19:24);
plot(smooth(offdiag,'sgolay'))
title('Off Diagonal Correlation');


end

