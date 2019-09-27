% sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/*.mat');
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat'); 
filter = 'WT';     
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
allCellsDES = nan(count,5); % [drugFRdiff,cntrlFRdiff, drugFREffectScore, drugCorrEffectScore, spike_depth(k)]
allCellsCorrMatrix = nan(count,300,300);
allCellsCorrScoreCurve = nan(count,300);
%%
sampleRate = 50; %hz
min = 15;
secInMin = 60;
minBeforeInjx = 5;
timeFRsamples = min * secInMin * sampleRate;
pre_timeFRsamples = minBeforeInjx * secInMin * sampleRate;
post_timeFRsamples = timeFRsamples - pre_timeFRsamples;

allCellsTimeFRcircaControlInjx = nan(count,timeFRsamples);
allCellsTimeFRcircaKetamineInjx = nan(count,timeFRsamples);
halfhourSampleNum = 30 * secInMin * sampleRate;
allCellsTimeFR30minAfterKetamineInjx = nan(count,halfhourSampleNum+1);
z = 0;
%
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
            'trial','all_cellCorrScore','trials_corrTemplate', 'avg_all_cellCorrScore', 'avg_cell_fr','trial_ds', 'all_frTime');

        spatialIndx = ismember(cells_to_plot,allSpatialIndx.(seshStr));
       
       if ~isempty(all_fr)
           all_fr = all_fr(spatialIndx,:,:);
           all_frTime = all_frTime(spatialIndx,:,:);
           nCells = size(all_fr,1);
           
           allCellsFR(z+1:z+nCells,:,:) = all_fr(1:nCells,1:300,1:200);
           allCellsDES(z+1:z+nCells,:,:) = all_drugEffectScores(1:nCells,1:5);
           allCellsCorrMatrix(z+1:z+nCells,:,:) = all_corrmatrix(1:nCells,1:300,1:300);
           allCellsCorrScoreCurve(z+1:z+nCells,:,:) = all_cellCorrScore(1:nCells,1:300);
           
           first50indx = find(trial_ds==50, 1,'first');
           first100indx = find(trial_ds==100, 1,'first');
           allCellsTimeFRcircaControlInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first50indx-pre_timeFRsamples+1:first50indx+post_timeFRsamples);
           allCellsTimeFRcircaKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:first100indx+post_timeFRsamples);
           allCellsTimeFR30minAfterKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx:first100indx+halfhourSampleNum);
           
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

%% PLOT Histfitket effect on FR score drug vs control
figure(3);
clf; hold on;

threshold = 10;

drugFRdiff = allCellsDES(allCellsDES(:,1) < threshold & allCellsDES(:,1) > -threshold,1);
cntrlFRdiff = allCellsDES(allCellsDES(:,2) < threshold & allCellsDES(:,2) > -threshold,2);

histogram(abs(cntrlFRdiff),50,'FaceAlpha',1);
histogram(abs(drugFRdiff),50,'FaceAlpha',.2);

set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('FR score for Control vs Ketamine(%s)',filter))
xlabel('FR Effect Score:Ratio of FR change between Drug & Control')
ylabel('Number of Cells')
vline(0,'r')



%% PLOT Histfitket effect on FR score
figure(4);
clf;

threshold = 15;
ketFRscore = allCellsDES(allCellsDES(:,3) < threshold & allCellsDES(:,3) > -threshold,3);
histfit(ketFRscore,120)
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Distribution of Ketamine FR Effect Scores(%s)',filter))
xlabel('FR Effect Score:Ratio of FR change between Drug & Control')
ylabel('Number of Cells')
vline(0,'r')

%% Plot Distribution of Ketamine Correlation Scores
figure(5);
hold on;
histfit(allCellsDES((allCellsDES(:,4) < 2),4))
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


% %% PLOT SCATTER of Corr Scorr vs Cell Depth
% figure(6);
% hold on;
% %         scatter(allCellsDES(:,3),allCellsDES(:,2), 150,'filled','r');
% scatter(allCellsDES(allCellsDES(:,2) < 2,3),allCellsDES(allCellsDES(:,2) < 2,2), 150,'filled','r');
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('Ketamine Correlation Effect Score vs Cell Depth')
% xlabel('Distance from Tip of Probe')
% ylabel('Pre vs Post Ketamine Correlation Effect Score')
% 
% 
%     
% %% PLOT SCATTER of FR score vs Cell Depth
% figure(8);
% hold on;
% scatter(allCellsDES(:,3),allCellsDES(:,1), 150,'filled','r');
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% set(gcf,'Position',[100 100 1000 1000])
% title('Ketamine FR Effect Score vs Cell Depth')
% xlabel('Distance from Tip of Probe')
% ylabel('Ketamine FR Effect Score(postDrugFR-preDrugFR)')


%% PLOT Scatter of Ketamine Corr Score vs ket effect on FR score
figure(9);
clf;

%         scatter(all_drugEffectScores(:,1),all_drugEffectScores(:,2), 150,'filled','r');
scatter(allCellsDES(abs(allCellsDES(:,3)) < 15,3),allCellsDES(abs(allCellsDES(:,3)) < 15,4), 150,'filled','r');
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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot each cell's CorrScore Curve
% all_cellCorrScore = nan(nCells, numel(trials_corrTemplate:size(allCellsFR,2)));
all_cellCorrScore = nan(nCells, numel(1:size(allCellsFR,2)));
for i = 1:size(allCellsFR,1)
   singleCellallTrialsFR = squeeze(allCellsFR(i,:,:));
   trials_corrTemplate = 50;
   [drugCorrEffectScore, cellCorrScore, corrTemplate] = calculateCorrScore(singleCellallTrialsFR, trials_corrTemplate);
   all_cellCorrScore(i,:) = cellCorrScore;
%    plot(cellCorrScore)
%    fprintf('Cell: %d, DrugCorrEffectScore = %.3f\n',i,drugCorrEffectScore);
%    
%    pause
end
fprintf('done')


close all;
figure(12);
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

%%
smMtx = smooth(allCellsTimeFRcircaControlInjx,sampleRate*10,'moving');
[~,sortIdx] = sort(allCellsTimeFRcircaKetamineInjx(:,pre_timeFRsamples)); % sort just the first column
sortedmat = allCellsTimeFRcircaKetamineInjx(sortIdx,:);   
imagesc(sortedmat,[0,30]); colorbar;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Firing Rate over Time 5 min before injection and 10 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
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

plot(smoothdata(nanmean(allCellsTimeFRcircaControlInjx,1),'gaussian',sampleRate*7),'LineWidth',2,'DisplayName','Control Injection');
plot(smoothdata(nanmean(allCellsTimeFRcircaKetamineInjx,1),'gaussian',sampleRate*7),'LineWidth',2,'DisplayName','Ketamine Injection');
% plot(smoothdata(nanmean(allCellsTimeFRcircaKetamineInjx,1),sampleRate*10,'moving'),'LineWidth',2,'DisplayName','Ketamine Injection');

% modify labels for tick marks
scaling  = sampleRate * secInMin; 
set(gca,'XLim',[0 15.1*scaling],'XTick',[0:scaling:15*scaling])
xticks = get(gca,'xtick');
x = 0:1:15;
newlabels = arrayfun(@(x) sprintf('%d', (x/scaling)-5), x*scaling, 'un', 0);
set(gca,'xticklabel',newlabels);
h = vline(5*scaling,'k','Injection');
title(sprintf('Control vs Ketamine Injections(%s)',filter))
xlabel('Minutes since Injection')
ylabel('Firing Rate (Hz)')
legend;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Firing Rate over Time 30 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
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

plot(smoothdata(nanmean(allCellsTimeFR30minAfterKetamineInjx,1),'gaussian',sampleRate*25),'LineWidth',2,'DisplayName','Ketamine Injection');
% plot(smoothdata(nanmean(allCellsTimeFRcircaKetamineInjx,1),sampleRate*10,'moving'),'LineWidth',2,'DisplayName','Ketamine Injection');
timeAfterDrug = 30;
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

%%
fs= 50;
pwelch(allCellsTimeFR60minAfterKetamineInjx(5,:), 10*fs, [] , [], fs)
pwelch(allCellsTimeFRcircaKetamineInjx(5,:), 10*fs, [] , [], fs)