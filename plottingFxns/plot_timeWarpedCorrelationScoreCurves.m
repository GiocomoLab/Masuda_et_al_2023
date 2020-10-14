function plot_timeWarpedCorrelationScoreCurves(cells)

% Load Params
if ~exist('paramsPath','var')
    params = readtable('./UniversalParams.xlsx');
else
    params = readtable(paramsPath);
end

samplingRate = params.TimeBin;
ds_factor = params.ds_factor;
secsInMin = 60; 
conversionFactor = secsInMin/(ds_factor*samplingRate);

%%
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
%% Normalized Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();  

avgControlInjxCorr = nanmean(all_cellCorrScore(:,1:50),2);
normalizedCorrScoreCurves = all_cellCorrScore./avgControlInjxCorr;

tw = timewarpTrialBasedScores(cells, normalizedCorrScoreCurves);

normCSC_data.y = convertMultidimensionalCellArray2paddedMatrix(tw.timewarpedScore)';
x = 1:size(normCSC_data.y,2);
normCSC_data.x = x./conversionFactor;
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
g.stat_summary('setylim','true');
g.set_names('x','Time','y','Correlation compared to Baseline Template (rho)','size',20); 
g.set_title('All','fontSize',30);
g.axe_property('FontSize',25);
g.set_color_options('chroma',0,'lightness',30);
g.draw;

% Baseline
figure();
baselineIndx_timewarpedScore = indexTimewarpedScoreOnBaselineOnline(tw,15);
normCSC_data.y = convertMultidimensionalCellArray2paddedMatrix(baselineIndx_timewarpedScore)';
normCSC_data.x = 1:size(normCSC_data.y,2);
normCSC_data.x = normCSC_data.x./conversionFactor;
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
g.stat_summary('setylim','true');
g.set_names('x','Time','y','Correlation compared to Baseline Template (rho)','size',20); 
g.set_title('Baseline','fontSize',30);
g.axe_property('FontSize',25);
g.set_color_options('chroma',0,'lightness',30);
g.draw;

% Control
figure();
cntrlIndx_timewarpedScore = indexTimewarpedScoreOnControlIndx(tw,15);
normCSC_data.y = convertMultidimensionalCellArray2paddedMatrix(cntrlIndx_timewarpedScore)';
normCSC_data.x = 1:size(normCSC_data.y,2);
normCSC_data.x = normCSC_data.x./conversionFactor;
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
g.stat_summary('setylim','true');
g.set_names('x','Time','y','Correlation compared to Baseline Template (rho)','size',20); 
g.set_title('Post-Control','fontSize',30);
g.axe_property('FontSize',25);
g.set_color_options('chroma',0,'lightness',30);
g.draw;

% Ketamine
figure();
ketIndx_timewarpedScore = indexTimewarpedScoreOnKetamineIndx(tw,15);
normCSC_data.y = convertMultidimensionalCellArray2paddedMatrix(ketIndx_timewarpedScore)';
normCSC_data.x = 1:size(normCSC_data.y,2);
normCSC_data.x = normCSC_data.x./conversionFactor;
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
g.stat_summary('setylim','true');
g.set_names('x','Time','y','Correlation compared to Baseline Template (rho)','size',20); 
g.set_title('Post-Ketamine','fontSize',30);
g.axe_property('FontSize',25);
g.set_color_options('chroma',0,'lightness',30);
g.draw;


end