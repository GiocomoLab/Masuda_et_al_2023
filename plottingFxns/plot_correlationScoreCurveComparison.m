function plot_correlationScoreCurveComparison(allCells, WTcellsIndx, KOcellsIndx, title)

WTCellsFR = allCells.spatialFR(WTcellsIndx,:,:);
KOCellsFR = allCells.spatialFR(KOcellsIndx,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot each cell's CorrScore Curve

trials_corrTemplate = 50;

WT_nCells = size(WTCellsFR,1);
WTall_cellCorrScore = nan(WT_nCells, numel(1:size(WTCellsFR,2)));
for i = 1:size(WTCellsFR,1)
   WTsingleCellallTrialsFR = squeeze(WTCellsFR(i,:,:));
   trials_corrTemplate = 50;
   [~, WTcellCorrScore, ~] = calculateCorrScore(WTsingleCellallTrialsFR, trials_corrTemplate);
   WTall_cellCorrScore(i,:) = WTcellCorrScore;
end

KO_nCells = size(KOCellsFR,1);
KOall_cellCorrScore = nan(KO_nCells, numel(1:size(KOCellsFR,2)));
for i = 1:size(KOCellsFR,1)
   KOsingleCellallTrialsFR = squeeze(KOCellsFR(i,:,:));

   [~, KOcellCorrScore, ~] = calculateCorrScore(KOsingleCellallTrialsFR, trials_corrTemplate);
   KOall_cellCorrScore(i,:) = KOcellCorrScore;
end
fprintf('done')

% 
% figure(); hold on;
% plot(nanmean(WTall_cellCorrScore,1),'LineWidth',5)
% plot(nanmean(KOall_cellCorrScore,1),'LineWidth',5)
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('Average Correlation Score Curve: WT vs KO')
% xlabel('Trial')
% ylabel('Correlation compared to Baseline Template (rho)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized Correlation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();  hold on;

WTavgControlInjxCorr = nanmean(nanmean(WTall_cellCorrScore(:,1:50),1));
KOavgControlInjxCorr = nanmean(nanmean(KOall_cellCorrScore(:,1:50),1));
avgControlInjxCorr = mean([WTavgControlInjxCorr,KOavgControlInjxCorr]);


WTnormalizedCorrScoreCurves = WTall_cellCorrScore/WTavgControlInjxCorr;
KOnormalizedCorrScoreCurves = KOall_cellCorrScore/KOavgControlInjxCorr;

plotTrialRange = 51:290;
normCSC_data.y = vertcat(WTnormalizedCorrScoreCurves(:,plotTrialRange),KOnormalizedCorrScoreCurves(:,plotTrialRange));
normCSC_data.x = repmat(plotTrialRange,size(normCSC_data.y,1),1);
normCSC_data.z = vertcat(repmat({'WT'},[WT_nCells,1]),repmat({'KO'},[KO_nCells,1]));

clf;
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y, 'color',normCSC_data.z);
g.stat_summary('type','sem','setylim','true');
g.set_names('x','Trial','y','Normalized Correlation to Baseline Template (rho)','size',20);
g.set_title(sprintf('Average Correlation Score Curve: %s',title),'FontSize',30);
g.axe_property('FontSize',25);
box off;
axis off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000]);
g.set_order_options('color',{'WT','KO'});
g.draw;


end