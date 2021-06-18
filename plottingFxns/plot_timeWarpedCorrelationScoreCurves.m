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
cellsFR = cells.spatialFRsmooth;
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

% avgControlInjxCorr = nanmean(all_cellCorrScore(:,1:50),2);
% normalizedCorrScoreCurves = all_cellCorrScore./avgControlInjxCorr;

tw = timewarpTrialBasedScores(cells, all_cellCorrScore);

normCSC_data.y = convertMultidimensionalCellArray2paddedMatrix(tw.timewarpedScore)';
x = 1:size(normCSC_data.y,2);
normCSC_data.x = x./conversionFactor;
normCSC_data.color = tw.session;

g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
% g.geom_line();
g.stat_summary('setylim','true');
g.set_names('x','Time','y','Correlation compared to Baseline Template (rho)','size',20); 
g.set_title('All','fontSize',30);
g.axe_property('FontSize',25);
% g.set_color_options('chroma',0,'lightness',30);
g.draw;

%% Ketamine
figure();
ketIndx_timewarpedScore = indexTimewarpedScoreOnKetamineIndx(tw,[]);
ket_y = convertMultidimensionalCellArray2paddedMatrix(ketIndx_timewarpedScore)';

normCSC_data.y = ket_y;
normCSC_data.x = 1:size(normCSC_data.y,2);
normCSC_data.x = normCSC_data.x./conversionFactor;

g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
g.stat_summary('setylim','true');
g.set_names('x','Time','y','Correlation compared to Baseline Template (rho)','size',15); 
g.set_title('Ketamine Correlation Score Curves','fontSize',30);
g.axe_property('FontSize',15);
g.set_color_options('chroma',0,'lightness',30);
g.draw;

%% Control
figure();
cntrlIndx_timewarpedScore = indexTimewarpedScoreOnControlIndx(tw,[]);
normCSC_data.y = convertMultidimensionalCellArray2paddedMatrix(cntrlIndx_timewarpedScore)';
normCSC_data.x = 1:size(normCSC_data.y,2);
normCSC_data.x = normCSC_data.x./conversionFactor;

g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
g.stat_summary('setylim','true');
g.set_names('x','Time','y','Correlation compared to Baseline Template (rho)','size',15); 
g.set_title('Control Correlation Score Curves','fontSize',30);
g.axe_property('FontSize',15);
g.set_color_options('chroma',0,'lightness',30);
g.draw;

%% Combined baseline, cntrl, ketamine
figure();
min = 10;
% baselineIndx_timewarpedScore = indexTimewarpedScoreOnBaselineOnline(tw,15);
% baseline_y = convertMultidimensionalCellArray2paddedMatrix(baselineIndx_timewarpedScore)';
cntrlIndx_timewarpedScore = indexTimewarpedScoreOnControlIndx(tw,min);
cntrl_y = convertMultidimensionalCellArray2paddedMatrix(cntrlIndx_timewarpedScore)';
ketIndx_timewarpedScore = indexTimewarpedScoreOnKetamineIndx(tw,min);
ket_y = convertMultidimensionalCellArray2paddedMatrix(ketIndx_timewarpedScore)';

normCSC_data.y = vertcat(cntrl_y,ket_y);
normCSC_data.x = 1:size(normCSC_data.y,2);
normCSC_data.x = normCSC_data.x./conversionFactor;
normCSC_data.color = vertcat(repmat({'Control'},size(cntrl_y,1),1),repmat({'Ketamine'},size(ket_y,1),1));


g=gramm('x',normCSC_data.x ,'y',normCSC_data.y,'color',normCSC_data.color);
g.stat_summary('setylim','true','type','ci');
g.set_names('x','Time','y','Correlation compared to Baseline Template (rho)','size',15); 
g.set_title('Correlation Score Curves','fontSize',30);
g.axe_property('FontSize',15);
% g.set_color_options('chroma',0,'lightness',30);
g.draw;


%% Stats
figure();

cntrlIndx_timewarpedScore = indexTimewarpedScoreOnControlIndx(tw,10);
cntrl_y = convertMultidimensionalCellArray2paddedMatrix(cntrlIndx_timewarpedScore)';
ketIndx_timewarpedScore = indexTimewarpedScoreOnKetamineIndx(tw,10);
ket_y = convertMultidimensionalCellArray2paddedMatrix(ketIndx_timewarpedScore)';

normCSC_data.y = vertcat(cntrl_y,ket_y);
normCSC_data.x = 1:size(normCSC_data.y,2);
normCSC_data.x = normCSC_data.x./conversionFactor;
normCSC_data.color = vertcat(repmat({'Control'},size(cntrl_y,1),1),repmat({'Ketamine'},size(ket_y,1),1));

[~,p] = ttest2(cntrl_y, ket_y);
clear g;
g=gramm('x',normCSC_data.color ,'y',normCSC_data.y,'color',normCSC_data.color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',0.15);
g.set_names('x','Time','y','Correlation compared to Baseline Template (rho)','size',15); 
g.set_title(sprintf('Corr Score p = %f',mean(p)),'fontSize',30);
g.axe_property('FontSize',15);
% g.set_color_options('chroma',0,'lightness',30);
g.draw;



end