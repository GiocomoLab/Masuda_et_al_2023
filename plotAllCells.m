%% Make Plots. Use after running Pool All cells
% Run all the plotting functions
% input: allCells struct

% add plotting functions to path
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('./plottingFxns'))

%% Generate Index 
WTcellsIndx = strcmp({allCells.metadata{:,4}}, 'WT')';
KOcellsIndx = strcmp({allCells.metadata{:,4}}, 'KO')';

stabilityTable = findStableCells(allCells); % {'totalStability', 'baselineStability', 'acuteDrugStability', 'endingStability', 'gainStability'}
stabilityThreshold = 0.45;
totalStabilityIndx = stabilityTable.totalStability > stabilityThreshold;
WTtotalStabilityIndx = (stabilityTable.totalStability > stabilityThreshold) & WTcellsIndx;
KOtotalStabilityIndx = (stabilityTable.totalStability > stabilityThreshold) & KOcellsIndx;

baselineStabilityIndx = stabilityTable.baselineStability > stabilityThreshold;
WTbaselineStabilityIndx = (stabilityTable.baselineStability > stabilityThreshold) & WTcellsIndx;
KObaselineStabilityIndx = (stabilityTable.baselineStability > stabilityThreshold) & KOcellsIndx;
%% Plot Stability Table STats
plot_stabilityTableStats(allCells, stabilityTable);

%% PLOT Histfit on FR score
plot_HistfitFRscore(allCells, WTcellsIndx,'WT');
plot_HistfitFRscore(allCells, KOcellsIndx,'KO');


plot_HistfitFRscore(allCells, WTtotalStabilityIndx,'WT');
plot_HistfitFRscore(allCells, KOtotalStabilityIndx,'KO');
%% Plot Distribution of Ketamine Correlation Scores
plot_HistfitKetCorrEffectScore(allCells, WTcellsIndx,'WT')
plot_HistfitKetCorrEffectScore(allCells, KOcellsIndx,'KO')

%% Plot Correlation Score Curves
plot_correlationScoreCurves(allCells, WTcellsIndx,'WT')
plot_correlationScoreCurves(allCells, KOcellsIndx,'KO')
plot_correlationScoreCurveComparison(allCells, WTcellsIndx, KOcellsIndx)

%% Plot Correlation Score Curve Comparisions
plot_correlationScoreCurveComparison(allCells, WTcellsIndx, KOcellsIndx, 'All Cells')
plot_correlationScoreCurveComparison(allCells, WTtotalStabilityIndx, KOtotalStabilityIndx, 'Whole Session Stable')
plot_correlationScoreCurveComparison(allCells, WTbaselineStabilityIndx, KObaselineStabilityIndx, 'Baseline Stable Cells')

%% Plot Correlation Score Curves by Animal

%% Plot Stability Score Curve
plot_stabilityScore(allCells, [], 'all cells')
plot_stabilityScore(allCells, WTcellsIndx,'WT')
plot_stabilityScore(allCells, KOcellsIndx,'KO')


%% Plot Correlation Matrix
plot_correlationMatrix(allCells,WTcellsIndx, 'WT')



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRELATION MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = fieldnames(allSpatialIndx);

sessionMap = [];
count = 0;
for k=1:numel(fn)
    if(isnumeric(allSpatialIndx.(fn{k})))
        seshCellNum = size(allSpatialIndx.(fn{k}),2);
        sessionMap(count+1:count+seshCellNum,1) = k;
        count = count + seshCellNum;
    end
end
%%
for k=1:numel(fn)
    try
    testCells = allCells.spatialFR(sessionMap == k,:,:);

    
    numCells = size(testCells,1);
    trialNum = size(testCells,2);
    spatialBins = size(testCells,3);

    flatFR2 = reshape(permute(testCells,[3 1 2]), [numCells*spatialBins, trialNum]);

    %P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
    corrMatrix = corr(fillmissing(flatFR2,'linear')); 
    figure(1); clf;
    imagesc(corrMatrix); colorbar;
    title(fn{k},'Interpreter', 'none');
    set(gca,'TickDir','out');
    set(gca,'ticklength',[0.005 0.025]);
    set(gca,'layer','bottom');
    box off;
    axis square;
    set(gca,'FontSize',30);
    set(gca,'FontName','Helvetica');
%     set(gcf,'Position',[100 100 1000 1000])

    tempAllCellScorrMatrix = nan(numCells, trialNum,trialNum);
    for i = 1:size(testCells,1)
        %calculate trial by trial correlation matrix for one cell 
        singleCellFR = squeeze(testCells(i,:,:))';
        corrMatrix = corr(fillmissing(singleCellFR,'linear'));
        tempAllCellScorrMatrix(i, :, :) = corrMatrix;
%         imagesc(corrMatrix)
%         pause 
    end
    figure(2); clf;
    imagesc(squeeze(nanmean(tempAllCellScorrMatrix,1))); colorbar;
    title(fn{k},'Interpreter', 'none');
    set(gca,'TickDir','out');
    set(gca,'ticklength',[0.005 0.025]);
    set(gca,'layer','bottom');
    box off;
    axis square;
    set(gca,'FontSize',30);
    set(gca,'FontName','Helvetica');
%     set(gcf,'Position',[100 100 1000 1000])
   pause 
    end
end

%%

testCells = allCells.spatialFR;


numCells = size(testCells,1);
trialNum = size(testCells,2);
spatialBins = size(testCells,3);

flatFR2 = reshape(permute(testCells,[3 1 2]), [numCells*spatialBins, trialNum]);

%P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
corrMatrix = corr(fillmissing(flatFR2,'linear')); 
figure(1); clf;
imagesc(corrMatrix); colorbar;
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box on;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Trial by Trial Population Activity Correlation Matrix(%s)',filter))


figure(2); clf;
imagesc(squeeze(nanmean(allCellsCorrMatrix,1)),[0, 0.2]); colorbar; 
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box on;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Avg Cell Trial by Trial Correlation Matrix(%s)',filter))


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Firing Rate over Time 5 min before injection and 10 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1); clf;
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.015 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% hold on;
% % plot(allCellsTimeFRcircaKetamineInjx)
% 
% plot(smoothdata(nanmean(allCellsTimeFRcircaControlInjx,1),'gaussian',sampleRate*7),'LineWidth',2,'DisplayName','Control Injection');
% plot(smoothdata(nanmean(allCellsTimeFRcircaKetamineInjx,1),'gaussian',sampleRate*7),'LineWidth',2,'DisplayName','Ketamine Injection');
% % plot(smoothdata(nanmean(allCellsTimeFRcircaKetamineInjx,1),sampleRate*10,'moving'),'LineWidth',2,'DisplayName','Ketamine Injection');
% 
% % modify labels for tick marks
% scaling  = sampleRate * secInMin; 
% set(gca,'XLim',[0 15.1*scaling],'XTick',[0:scaling:15*scaling])
% xticks = get(gca,'xtick');
% x = 0:1:15;
% newlabels = arrayfun(@(x) sprintf('%d', (x/scaling)-5), x*scaling, 'un', 0);
% set(gca,'xticklabel',newlabels);
% h = vline(5*scaling,'k','Injection');
% title(sprintf('Control vs Ketamine Injections(%s)',filter))
% xlabel('Minutes since Injection')
% ylabel('Firing Rate (Hz)')
% legend;

%
figure(); clf;
clear g;
xsize = size(allCellsTimeFRcircaControlInjx,2);
numCell = size(allCellsTimeFRcircaControlInjx,1);
dsfactor = 50;
x = downsample(1:xsize,dsfactor);
x = x/scaling - 5;
y = vertcat(allCellsTimeFRcircaControlInjx,allCellsTimeFRcircaKetamineInjx);
y = downsample(y',dsfactor)';
% z = vertcat(ones([numCell,1]),zeros(numCell,1)); %control vs ket label
z = cellstr(vertcat(repmat("Control",[numCell,1]),repmat("Ketamine",[numCell,1])));


g(1,1) = gramm('x',x,'y',y,'color',z); 
g(1,1).stat_summary('type','sem','setylim',true);

g(1,1).set_title(sprintf('Control vs Ketamine Injections(%s)',filter), 'FontSize', 40);
g(1,1).set_names('x','Minutes since Injection','y','Firing Rate (Hz)');
g(1,1).set_text_options('base_size',20);
g.draw()
axis square;
set(gcf,'Position',[100 100 1000 1000])
set(gca,'TickDir','out');
a

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Firing Rate over Time 45 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(); clf;
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
hold on;
% plot(allCellsTimeFRcircaKetamineInjx)

plot(smoothdata(nanmean(allCellsTimeFR45minAfterKetamineInjx,1),'gaussian',sampleRate*25),'LineWidth',2,'DisplayName','Ketamine Injection');
% plot(smoothdata(nanmean(allCellsTimeFRcircaKetamineInjx,1),sampleRate*10,'moving'),'LineWidth',2,'DisplayName','Ketamine Injection');
timeAfterDrug = 45;
% modify labels for tick marks
scaling  = sampleRate * secInMin;
tickSteps = 5;
set(gca,'XLim',[0 timeAfterDrug*scaling],'XTick',[0:tickSteps*scaling:timeAfterDrug*scaling])
xticks = get(gca,'xtick');
x = 0:tickSteps:timeAfterDrug;
newlabels = arrayfun(@(x) sprintf('%d', x/scaling), x*scaling, 'un', 0);
set(gca,'xticklabel',newlabels);
title(sprintf('Ketamine-induced FR over Time(%s)',filter))
xlabel('Minutes since Injection')
ylabel('Firing Rate (Hz)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Firing Rate over Time 5min before and 60 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(); clf;
clear g;
scaling  = sampleRate * secInMin;
downsampleFactor = 250;
timeAfterDrug = 60;
timeBeforeDrug = -5;
x = timeBeforeDrug*scaling+1:downsampleFactor:timeAfterDrug*scaling;
y = downsample(allCells.timeFRneg5to60minAfterKetamineInjx',downsampleFactor)';

x_min = x./scaling;0

g(1,1) = gramm('x',x_min,'y',y); 
g(1,1).stat_summary('type','sem');
% g(1,1).stat_smooth();

% g(1,1).geom_line();
g(1,1).axe_property('YLim',[5 10]);
g(1,1).set_title(sprintf('Ketamine-induced FR over Time(%s)',filter), 'FontSize', 40);
g(1,1).set_names('x','Minutes since Injection','y','Firing Rate(Hz)');
g(1,1).set_text_options('base_size',20);
g.draw()
axis square;
set(gcf,'Position',[100 100 1000 1000])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Firing Rate over Time 45 min after injection with KO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(); clf; hold on;
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
hold on;
% plot(allCellsTimeFRcircaKetamineInjx)

plot(smoothdata(nanmean(allCellsTimeFR45minAfterKetamineInjx,1),'gaussian',sampleRate*25),'LineWidth',2,'DisplayName','WT');
plot(smoothdata(nanmean(allCellsTimeFR45minAfterKetamineInjxKO,1),'gaussian',sampleRate*25),'LineWidth',2,'DisplayName','HCN1ko');
% plot(smoothdata(nanmean(allCellsTimeFRcircaKetamineInjx,1),sampleRate*10,'moving'),'LineWidth',2,'DisplayName','Ketamine Injection');
timeAfterDrug = 45;
% modify labels for tick marks
scaling  = sampleRate * secInMin;
tickSteps = 5;
set(gca,'XLim',[0 timeAfterDrug*scaling],'XTick',[0:tickSteps*scaling:timeAfterDrug*scaling])
xticks = get(gca,'xtick');
x = 0:tickSteps:timeAfterDrug;
newlabels = arrayfun(@(x) sprintf('%d', x/scaling), x*scaling, 'un', 0);
set(gca,'xticklabel',newlabels);
title(sprintf('Ketamine-induced FR over Time(WT vs HCN1ko)',filter))
xlabel('Minutes since Injection')
ylabel('Firing Rate (Hz)')
legend;

