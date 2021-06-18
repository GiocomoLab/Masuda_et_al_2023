function normCSC_data = plot_timeWarpedStabilityScoreCurves(cells)

addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('/Users/keimasuda/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
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
all_cellStabilityScore = nan(nCells, numel(1:size(cellsFR,2)));
for i = 1:size(cellsFR,1)
   singleCellallTrialsFR = squeeze(cellsFR(i,:,:));
   stabilityScore = calculateStabilityScore(singleCellallTrialsFR);
   all_cellStabilityScore(i,:) = stabilityScore;
end
fprintf('done\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalized Stability Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
tw = timewarpTrialBasedScores(cells, all_cellStabilityScore);

% 
figure();  


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
g.set_color_options('chroma',0,'lightness',30);
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
g.set_names('x','Time','y','Stability(rho)','size',15); 
g.set_title('Ketamine Stability Score Curves','fontSize',30);
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
g.set_names('x','Time','y','Stability (rho)','size',15); 
g.set_title('Control Stability Score Curves','fontSize',30);
g.axe_property('FontSize',15);
g.set_color_options('chroma',0,'lightness',30);
g.draw;

%% Combined cntrl, ketamine 10 minutes
close all;

min = 10;
baselineIndx_timewarpedScore = indexTimewarpedScoreOnBaselineOnline(tw,min);
baseline_y = convertMultidimensionalCellArray2paddedMatrix(baselineIndx_timewarpedScore)';
cntrlIndx_timewarpedScore = indexTimewarpedScoreOnControlIndx(tw,min);
cntrl_y = convertMultidimensionalCellArray2paddedMatrix(cntrlIndx_timewarpedScore)';
ketIndx_timewarpedScore = indexTimewarpedScoreOnKetamineIndx(tw,min);
ket_y = convertMultidimensionalCellArray2paddedMatrix(ketIndx_timewarpedScore)';

normCSC_data.y = vertcat(baseline_y,cntrl_y,ket_y);
normCSC_data.x = 1:size(normCSC_data.y,2);
normCSC_data.x = normCSC_data.x./conversionFactor;
normCSC_data.color = vertcat(repmat({'Baseline'},size(baseline_y,1),1),repmat({'Control'},size(cntrl_y,1),1),repmat({'Ketamine'},size(ket_y,1),1));


figure();
clear g;
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y,'color',normCSC_data.color);
g.stat_summary('setylim','true','type','bootci');
g.set_color_options('chroma',40)
g.set_names('x','Time','y','Stability (rho)','size',15); 
g.set_title('Stability Score Curves','fontSize',30);
g.axe_property('FontSize',15);
% g.set_color_options('chroma',0,'lightness',30);
g.draw;


%% Combined cntrl, ketamine 60 minutes
close all;

min = 30;
baselineIndx_timewarpedScore = indexTimewarpedScoreOnBaselineOnline(tw,min);
baseline_y = convertMultidimensionalCellArray2paddedMatrix(baselineIndx_timewarpedScore)';
cntrlIndx_timewarpedScore = indexTimewarpedScoreOnControlIndx(tw,min);
cntrl_y = convertMultidimensionalCellArray2paddedMatrix(cntrlIndx_timewarpedScore)';
ketIndx_timewarpedScore = indexTimewarpedScoreOnKetamineIndx(tw,min);
ket_y = convertMultidimensionalCellArray2paddedMatrix(ketIndx_timewarpedScore)';

normCSC_data.y = vertcat(baseline_y,cntrl_y,ket_y);
normCSC_data.x = 1:size(normCSC_data.y,2);
normCSC_data.x = normCSC_data.x./conversionFactor;
normCSC_data.color = vertcat(repmat({'Baseline'},size(baseline_y,1),1),repmat({'Control'},size(cntrl_y,1),1),repmat({'Ketamine'},size(ket_y,1),1));


figure();
clear g;
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y,'color',normCSC_data.color);
g.stat_summary('setylim','true','type','bootci');
g.set_color_options('chroma',40)
g.set_names('x','Time','y','Stability (rho)','size',15); 
g.set_title('Stability Score Curves','fontSize',30);
g.axe_property('FontSize',15);
% g.set_color_options('chroma',0,'lightness',30);
g.draw;

%%

figure();
clear g;
g=gramm('x',normCSC_data.color ,'y',normCSC_data.y,'color',normCSC_data.color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',0.15);
g.set_names('x',[],'y','Stability (rho)','size',20); 
g.set_title(sprintf('Stability Score'),'fontSize',20);
g.axe_property('FontSize',12);
% g.set_color_options('chroma',0,'lightness',30);
g.draw;

% Stats
[~,p1] = ttest2(baseline_y, cntrl_y);
[~,p2] = ttest2(cntrl_y, ket_y);
[~,p3] = ttest2(baseline_y, ket_y);
fprintf('Baseline vs cntrl: p%0.03f\nCntrl vs Ket: p%0.03f \nBaseline vs ket: p%0.03f\n',mean(p1),mean(p2),mean(p3))

[p1,h] = ranksum(nanmean(baseline_y,2), nanmean(cntrl_y,2));
[p2,h] = ranksum(nanmean(cntrl_y,2), nanmean(ket_y,2));
[p3,h] = ranksum(nanmean(baseline_y,2), nanmean(ket_y,2));
fprintf('Ranksum: \nBaseline vs cntrl: p%0.03f\nCntrl vs Ket: p%0.03f \nBaseline vs ket: p%0.03f\n',p1,p2,p3)



end