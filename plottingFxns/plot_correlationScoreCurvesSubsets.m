function plot_correlationScoreCurvesSubsets(allCells,cellIndx)

% cellIndx = baselineStabilityIndx;
allCellsFR = allCells.spatialFR(cellIndx,:,:);
metadata = allCells.metadata(cellIndx,:);

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized  Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();


trials = 51:290;
avgControlInjxCorr = nanmean(nanmean(all_cellCorrScore(:,1:50),1));
normalizedCorrScoreCurves = all_cellCorrScore./avgControlInjxCorr;
normCSC_data.y = normalizedCorrScoreCurves(:,trials);
normCSC_data.x = repmat(trials,size(normalizedCorrScoreCurves,1),1);
excludeAnimal = metadata(:,2);
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y,'color',categorical(metadata(:,4)),'subset',~strcmp(excludeAnimal,['npI1']));

g.stat_summary('setylim','true');
g.set_names('x','Trial','y','Correlation compared to Baseline Template (rho)','size',20); 
g.set_title('Normalized Corr Scorr Curve by Animal','fontSize',30);
% g.axe_property('FontSize',25);
box off;
axis off;
% axis square;
% 
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
g.set_order_options('color',{'WT','KO'});
g.facet_wrap(categorical(metadata(:,2)),'scale','independent');
g.draw;


end