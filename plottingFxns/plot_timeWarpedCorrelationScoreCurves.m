function plot_timeWarpedCorrelationScoreCurves(cells, filter)

cellsFR = cells.spatialFR2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot each cell's CorrScore Curve
nCells = size(cellsFR,1);
all_cellCorrScore = nan(nCells, numel(1:size(cellsFR,2)));
for i = 1:size(cellsFR,1)
   singleCellallTrialsFR = squeeze(cellsFR(i,:,:));
   trials_corrTemplate = 50;
   [~, cellCorrScore, ~] = calculateCorrScore(singleCellallTrialsFR, trials_corrTemplate);
   all_cellCorrScore(i,:) = cellCorrScore;
end
fprintf('done\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalized  Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();  

avgControlInjxCorr = nanmean(all_cellCorrScore(:,1:50),2);
normalizedCorrScoreCurves = all_cellCorrScore./avgControlInjxCorr;

ds_factor = 100;
timewarpedScore = timewarpTrialBasedScores(cells, normalizedCorrScoreCurves,ds_factor);


normCSC_data.y = timewarpedScore;
normCSC_data.x = 1:size(timewarpedScore,2);
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
g.stat_summary('setylim','true');
g.set_names('x','Trial','y','Correlation compared to Baseline Template (rho)','size',20); 
g.set_title(sprintf('Average Correlation Score Curve(%s - %d Cells)',filter,size(cellsFR,1)),'fontSize',30);
g.axe_property('FontSize',25);
g.set_color_options('chroma',0,'lightness',30);
g.draw;


end