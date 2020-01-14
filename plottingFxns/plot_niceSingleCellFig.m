function plot_niceSingleCellFig(allCells)

%%
i = randi(size(allCells.spatialFR,1));
% i = 1148
% i = 300
% i = 3213
%
singleCellFR = squeeze(allCells.spatialFR(i,:,:));

lickx = allCells.lickX(i).lickx;
lickt = allCells.lickT(i).lickt;
posx = allCells.posX(i).posx;
post = allCells.posT(i).post;
FRtime = allCells.FRtime(i).FRtime';
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

linearFractionalOccupancy = calculate_1D_LFO(posx,post);
[I_sec, I_spike] = calculate_1DspatialInformation(singleCellFR,linearFractionalOccupancy);
sum(sum(singleCellFR<0))

trialsPerBlock = 10;
trialBlockSpatialInformation = nan(max(trial)/trialsPerBlock,2);
for j = 1:max(trial)/trialsPerBlock
    
    trialStart = (j-1)*trialsPerBlock+1;
    trialEnd = trialStart+trialsPerBlock-1;
    trialRange = trialStart:trialEnd;   
    
    trial_posx = posx( trial>=trialStart & trial<=trialEnd);
    trial_post = post( trial>=trialStart & trial<=trialEnd);
    linearFractionalOccupancy = calculate_1D_LFO(trial_posx,trial_post);
    
    [I_sec, I_spike] = calculate_1DspatialInformation(singleCellFR(trialRange,:),linearFractionalOccupancy);
    trialBlockSpatialInformation(j,:) = [I_sec, I_spike];
end
trialBlockSpatialInformation(isnan(trialBlockSpatialInformation))=0;

%
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

plot(trialBlockSpatialInformation(:,1))
scaling = trialsPerBlock;
newlabels = arrayfun(@(x) sprintf('%d', scaling * x), xticks, 'un', 0);
set(gca,'xticklabel',newlabels)

set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
% set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Spatial Information = Isec:%0.3f; Ispike:%0.3f', I_sec, I_spike))
% xlabel('Trial Number')
ylabel('SI Score')
% legend('I sec', 'I spike');

subplot(row,col,20:24);
% bar(lickAccuracyByTrial);
scatter(trial(lick_idx),posx(lick_idx));
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
% set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
title('Licks')
xlabel('Trial Number')
ylabel('VR cm')
end