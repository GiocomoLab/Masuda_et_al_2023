function plot_correlationScoreCurves(allCells, cellIndx, filter)

allCellsFR = allCells.spatialFR(cellIndx,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot each cell's CorrScore Curve
% all_cellCorrScore = nan(nCells, numel(trials_corrTemplate:size(allCellsFR,2)));
nCells = size(allCellsFR,1);
all_cellCorrScore = nan(nCells, numel(1:size(allCellsFR,2)));
for i = 1:size(allCellsFR,1)
   singleCellallTrialsFR = squeeze(allCellsFR(i,:,:));
   trials_corrTemplate = 50;
   [~, cellCorrScore, ~] = calculateCorrScore(singleCellallTrialsFR, trials_corrTemplate);
   all_cellCorrScore(i,:) = cellCorrScore;
%    plot(cellCorrScore)
%    fprintf('Cell: %d, DrugCorrEffectScore = %.3f\n',i,drugCorrEffectScore);
%    
%    pause
end
fprintf('done')


figure(1); hold on;
plot(nanmean(all_cellCorrScore,1),'LineWidth',5)
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('%d Cell, Average Correlation Score Curve(%s)',size(allCellsFR,1),filter))
xlabel('Trial')
ylabel('Corrleation compared to Baseline Template (rho)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized  Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);  hold on;

avgControlInjxCorr = nanmean(nanmean(all_cellCorrScore(:,1:50),1));
normalizedCorrScoreCurves = all_cellCorrScore./avgControlInjxCorr;
normCSC_data.y = normalizedCorrScoreCurves(:,51:300);
normCSC_data.x = repmat(51:300,size(normalizedCorrScoreCurves,1),1);
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
% g.geom_line();
g.stat_summary('setylim','true');
g.set_names('x','Trial','y','Correlation compared to Baseline Template (rho)','size',20); 
g.set_title(sprintf('Average Correlation Score Curve(%s - %d Cells)',filter,size(allCellsFR,1)),'fontSize',30);
g.axe_property('FontSize',25);
box off;
axis off;
axis square;

% normalizedCorrScoreCurve = nanmean(all_cellCorrScore,1)./avgControlInjxCorr;
% plot(51:300,normalizedCorrScoreCurve(51:end),'LineWidth',5','DisplayName','WT')
% % plot(51:300,normalizedCorrScoreCurveKO(51:end),'LineWidth',5,'DisplayName','HCN1ko')
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
% title(sprintf('Average Correlation Score Curve'))
% xlabel('Trial')
% ylabel('Corrleation Score')
% legend;

g.draw;


end