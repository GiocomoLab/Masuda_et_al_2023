% sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/*.mat');
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat'); 
filter = 'KO';     
sessions = filterSessions(sessions, filter);
load(sprintf('/Users/KeiMasuda/Desktop/fkm_analysis/allSpatialIndx%s_01.mat',filter));
% load(sprintf('/Users/KeiMasuda/Desktop/fkm_analysis/allIndx%s.mat',filter));

fn = fieldnames(allSpatialIndx);
count = 0;
for k=1:numel(fn)
    if(isnumeric(allSpatialIndx.(fn{k})))
        count = count + size(allSpatialIndx.(fn{k}),2);
    end
end

allCellsFR = nan(count,300,200);
allCellsDES = nan(count,3);
allCellsCorrMatrix = nan(count,300,300);
allCellsCorrScoreCurve = nan(count,251);
z = 0;

for n = 1:numel(fn)
    try
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        session_name = sessions(n).name(1:end-4);
        animalName = extractBefore(session_name,'_');
        sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
        seshStr = sprintf('%s_%s',animalName, sessionDate);
        trackLength = 400;
        load(fullfile(matPath), 'all_fr', 'avg_all_fr', 'all_corrmatrix', 'avg_all_corrmatrix', ...
             'all_waveforms', 'cells_to_plot','spike_depth','all_drugEffectScores',...
            'trial','all_cellCorrScore','trials_corrTemplate', 'avg_all_cellCorrScore', 'avg_cell_fr');

       spatialIndx = ismember(cells_to_plot,allSpatialIndx.(seshStr));
       
       if ~isempty(all_fr)
           all_fr = all_fr(spatialIndx,:,:);
           nCells = size(all_fr,1);
           
           allCellsFR(z+1:z+nCells,:,:) = all_fr(1:nCells,1:300,1:200);
           allCellsDES(z+1:z+nCells,:,:) = all_drugEffectScores(1:nCells,1:3);
           allCellsCorrMatrix(z+1:z+nCells,:,:) = all_corrmatrix(1:nCells,1:300,1:300);
           allCellsCorrScoreCurve(z+1:z+nCells,:,:) = all_cellCorrScore(1:nCells,1:251);
           
           z = z + nCells;
           fprintf('Session: %d; Adding %d for %d/%d cells\n', n,nCells,z,count)
       else
           fprintf('Adding 0 for %d/%d cells\n',z,count)
       end
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end
fprintf('done\n')
%%
plot(squeeze(nanmean(squeeze(nanmean(allCellsFR,1)),2)))


%% Plot Distribution of Ketamine Correlation Scores
figure(5);
hold on;
histfit(allCellsDES((allCellsDES(:,2) < 2),2))
%         histfit(allCellsDES(:,2),50, 'kernel')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Distribution of Ketamine Correlation Effect Scores(%s)',filter))
xlabel('Pre vs Post Ketamine Correlation Effect Score')
ylabel('Number of Cells')


%% PLOT SCATTER of Corr Scorr vs Cell Depth
figure(6);
hold on;
%         scatter(allCellsDES(:,3),allCellsDES(:,2), 150,'filled','r');
scatter(allCellsDES(allCellsDES(:,2) < 2,3),allCellsDES(allCellsDES(:,2) < 2,2), 150,'filled','r');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title('Ketamine Correlation Effect Score vs Cell Depth')
xlabel('Distance from Tip of Probe')
ylabel('Pre vs Post Ketamine Correlation Effect Score')


    
%% PLOT SCATTER of FR score vs Cell Depth
figure(8);
hold on;
scatter(allCellsDES(:,3),allCellsDES(:,1), 150,'filled','r');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title('Ketamine FR Effect Score vs Cell Depth')
xlabel('Distance from Tip of Probe')
ylabel('Ketamine FR Effect Score(postDrugFR-preDrugFR)')


%% PLOT Scatter of Ketamine Corr Score vs ket effect on FR score
figure(9);
clf;

%         scatter(all_drugEffectScores(:,1),all_drugEffectScores(:,2), 150,'filled','r');
scatter(allCellsDES(allCellsDES(:,2) < 2,1),allCellsDES(allCellsDES(:,2) < 2,2), 150,'filled','r');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Ketamine Correlation Effect Score vs Firing Rate Score(%s)',filter))
xlabel('FR Score Cell (Hz)')
ylabel('Pre vs Post Ketamine Correlation Effect Score')

%% PLOT Histfitket effect on FR score
figure(10);
clf;
%         boxplot(allCellsDES(:,1))
histfit(allCellsDES(:,1),50)
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Distribution of Ketamine FR Effect Scores(%s)',filter))
xlabel('Pre vs Post Ketamine FR Effect Score')
ylabel('Number of Cells')

%%
close all;
figure(11);
imagesc(squeeze(nanmean(allCellsCorrMatrix,1)),[0 0.15])
colorbar;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
title(sprintf('Trial by Trial Population Correlation Matrix(%s)',filter))


numCells = size(allCellsFR,1);
trialNum = size(allCellsFR,2);
spatialBins = size(allCellsFR,3);

flatFR2 = reshape(permute(allCellsFR,[3 1 2]), [numCells*spatialBins, trialNum]);
%P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
corrMatrix = corr(fillmissing(flatFR2,'linear')); 
figure(12);
imagesc(corrMatrix); colorbar;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
title(sprintf('Avg Cell Trial by Trial Correlation Matrix(%s)',filter))

%% Plot each cell's CorrScore Curve
% all_cellCorrScore = nan(nCells, numel(trials_corrTemplate:size(allCellsFR,2)));
all_cellCorrScore = nan(nCells, numel(1:size(allCellsFR,2)));
for i = 1:size(allCellsFR,1)
   singleCellallTrialsFR = squeeze(allCellsFR(i,:,:));
   trials_corrTemplate = 50;
   [drugCorrEffectScore, cellCorrScore, corrTemplate] = calculateCorrScore(singleCellallTrialsFR, trials_corrTemplate);
   all_cellCorrScore(i,:) = cellCorrScore;
   plot(cellCorrScore)
%    fprintf('Cell: %d, DrugCorrEffectScore = %.3f\n',i,drugCorrEffectScore);
%    
%    pause
end
fprintf('done')
%
close all;
figure(12);
plot(mean(all_cellCorrScore,1),'LineWidth',5)
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
    testCells = allCellsFR(sessionMap == k,:,:);

    
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

testCells = allCellsFR;


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