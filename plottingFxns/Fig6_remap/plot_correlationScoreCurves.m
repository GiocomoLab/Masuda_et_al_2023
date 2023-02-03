function plot_correlationScoreCurves(cells, filter)

cellsFR = cells.spatialFRsmooth;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot each cell's CorrScore Curve
% all_cellCorrScore = nan(nCells, numel(trials_corrTemplate:size(cellsFR,2)));
nCells = size(cellsFR,1);
all_cellCorrScore = nan(nCells, numel(1:size(cellsFR,2)));
for i = 1:size(cellsFR,1)
   singleCellallTrialsFR = squeeze(cellsFR(i,:,:));
   trials_corrTemplate = 50;
   [~, cellCorrScore, ~] = calculateCorrScore(singleCellallTrialsFR, trials_corrTemplate);
   all_cellCorrScore(i,:) = cellCorrScore;
%    plot(cellCorrScore)
%    fprintf('Cell: %d, DrugCorrEffectScore = %.3f\n',i,drugCorrEffectScore);
%    
%    pause
end
fprintf('done\n')

% 
% figure(1); hold on;
% plot_lineWithSEM(all_cellCorrScore,[])
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title(sprintf('%d Cell, Average Correlation Score Curve(%s)',size(cellsFR,1),filter))
% xlabel('Trial')
% ylabel('Corrleation compared to Baseline Template (rho)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized  Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();  hold on;
plotTrialRange = 51:290;
avgControlInjxCorr = nanmean(all_cellCorrScore(:,1:50),2);
normalizedCorrScoreCurves = all_cellCorrScore./avgControlInjxCorr;




normCSC_data.y = normalizedCorrScoreCurves(:,plotTrialRange);
normCSC_data.x = repmat(plotTrialRange,size(normalizedCorrScoreCurves,1),1);
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
% g.geom_line();
g.stat_summary('setylim','true');
g.set_names('x','Trial','y','Correlation compared to Baseline Template (rho)','size',20); 
g.set_title(sprintf('Average Correlation Score Curve(%s - %d Cells)',filter,size(cellsFR,1)),'fontSize',30);
g.axe_property('FontSize',25);
g.set_color_options('chroma',0,'lightness',30);
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