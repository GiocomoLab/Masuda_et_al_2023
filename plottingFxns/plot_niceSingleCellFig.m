function plot_niceSingleCellFig(allCells)

%%
i = randi(size(allCells.spatialFR,1));
% i = 1148
% i = 300
%%
singleCellFR = squeeze(allCells.spatialFR(i,:,:));

figure(1);
row = 4;
col = 6;

subplot(row,col,[1,7,13,19])
imagesc(singleCellFR);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
% set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Single Cell(%d)',i))
xlabel('VR cm')
ylabel('Trial Number')

%
speed = allCells.speed(i).speed;
trial = allCells.trial(i).trial;


meanSpeed = nan(max(trial),1);
for j = 1:max(trial)
    meanSpeed(j) = nanmean(speed(trial == j));
end
subplot(row,col,2:6)
plot(smooth(meanSpeed,'sgolay'));
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
% set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title('Speed')
ylabel('cm/s')
%
subplot(row,col,8:12)
stabilityScoreCurve = squeeze(allCells.stabilityScoreCurve(i,:));
plot(smooth(stabilityScoreCurve,'sgolay'));
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
% set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title('Stability Score')
ylabel('Stability Score')


subplot(row,col,14:18);
lickx = allCells.lickX(i).lickx;
lickt = allCells.lickT(i).lickt;
posx = allCells.posX(i).posx;
post = allCells.posT(i).post;
FRtime = allCells.FRtime(i).FRtime';
[~,~,lick_idx] = histcounts(lickt,post);

lickAccuracyByTrial = zeros(1,max(trial));
for i = 1:max(trial)
    trialLicks = lickx(trial(lick_idx) == i);
    goodLicks = sum(trialLicks<25) + sum(trialLicks>max(posx)-25); 
    if trialLicks ~= 0
        lickAccuracyByTrial(i) = goodLicks/numel(trialLicks);
    else
        lickAccuracyByTrial(i) = 0.0;
    end
end
plot(lickAccuracyByTrial);

set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
% set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title('lickAccuracy Score')
xlabel('Trial Number')
ylabel('VR cm')

subplot(row,col,20:24);
scatter(trial(lick_idx),posx(lick_idx));
end