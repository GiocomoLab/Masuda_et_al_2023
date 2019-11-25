% sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/*.mat');
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat'); 
filter = 'WT';     
sessions = filterSessions(sessions, filter);
% load(sprintf('/Users/KeiMasuda/Desktop/fkm_analysis/allSpatialIndx%s_01.mat',filter));
load(sprintf('/Users/KeiMasuda/Desktop/fkm_analysis/allSpatialIndx%s.mat',filter));

fn = fieldnames(allSpatialIndx);
count = 0;
seshVector = [];
for k=1:numel(fn)
    if(isnumeric(allSpatialIndx.(fn{k})))
        count = count + size(allSpatialIndx.(fn{k}),2);
        seshVector = vertcat(seshVector, repmat(k,size(allSpatialIndx.(fn{k}),2),1));
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
fortyfiveminSampleNum = 45 * secInMin * sampleRate;
fiftyMinSampleNum = 50 * secInMin * sampleRate;
sixtyMinSampleNum = 60 * secInMin * sampleRate;
sixtyfiveMinSampleNum = 65 * secInMin * sampleRate;
allCellsTimeFR30minAfterKetamineInjx = nan(count,halfhourSampleNum+1);
allCellsTimeFR45minAfterKetamineInjx = nan(count,fortyfiveminSampleNum+1);
allCellsTimeFRneg5to45minAfterKetamineInjx = nan(count,fiftyMinSampleNum);
allCellsTimeFRneg5to60minAfterKetamineInjx = nan(count,sixtyfiveMinSampleNum);
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
%            allCellsTimeFR30minAfterKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx:first100indx+halfhourSampleNum);
%            if first100indx+fortyfiveminSampleNum < size(all_frTime,2)
%                allCellsTimeFR45minAfterKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx:first100indx+fortyfiveminSampleNum);
%                allCellsTimeFRneg5to45minAfterKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:first100indx+fortyfiveminSampleNum);
%            else
%                time45IndexLength = 1:size(all_frTime(1:nCells,first100indx:end),2);
%                maxTimeIndexLength = 1:size(all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:end),2);
%                
%                cellIndexLength = z+1:z+nCells;
%                allCellsTimeFR45minAfterKetamineInjx(cellIndexLength,timeIndexLength) = all_frTime(1:nCells,first100indx:end);
%                allCellsTimeFRneg5to45minAfterKetamineInjx(cellIndexLength,maxTimeIndexLength) = all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:end);
%                
%            end
           
           if first100indx+sixtyMinSampleNum < size(all_frTime,2)
               allCellsTimeFRneg5to60minAfterKetamineInjx(z+1:z+nCells,:) = all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:first100indx+sixtyMinSampleNum);
           else
               maxTimeIndexLength = 1:size(all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:end),2);
               cellIndexLength = z+1:z+nCells;
               allCellsTimeFRneg5to60minAfterKetamineInjx(cellIndexLength,maxTimeIndexLength) = all_frTime(1:nCells,first100indx-pre_timeFRsamples+1:end); 
           end
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



%% PLOT Histfit on FR score
figure(4);
clf;

threshold = 15;
ketFRscore = allCellsDES(allCellsDES(:,3) < threshold & allCellsDES(:,3) > -threshold,3);
histfit(ketFRscore,50)
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
histfit(allCellsDES((allCellsDES(:,4) < 2),4),50)
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
vline(0,'r')

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
xlabel('FR Score Cell')
ylabel('Pre vs Post Ketamine Correlation Effect Score')


%%
close all;
figure(11);
imagesc(squeeze(nanmean(allCellsCorrMatrix,1)),[0 0.3])
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
% Noramlized  Correleation Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(13);

avgControlInjxCorr = nanmean(nanmean(all_cellCorrScore(:,1:50),1));
normalizedCorrScoreCurves = all_cellCorrScore./avgControlInjxCorr;
normCSC_data.y = normalizedCorrScoreCurves(:,51:300);
normCSC_data.x = repmat(51:300,size(normalizedCorrScoreCurves,1),1);
g=gramm('x',normCSC_data.x ,'y',normCSC_data.y);
% g.geom_line();
g.stat_summary('setylim','true');
g.set_names('x','Trial','y','Correlation compared to Baseline Template (rho)','size',20); 
g.set_title(sprintf('Average Correlation Score Curve(%s - %d Cells)',filter,size(allCellsFR,1)),'fontSize',30);
g.axe_property('FontSize',25)
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
% Plot Firing Rate over Time 45 min after injection
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

plot(smoothdata(nanmean(allCellsTimeFR45minAfterKetamineInjx,1),'gaussian',sampleRate*25),'LineWidth',2,'DisplayName','Ketamine Injection');
% plot(smoothdata(nanmean(allCellsTimeFRcircaKetamineInjx,1),sampleRate*10,'moving'),'LineWidth',2,'DisplayName','Ketamine Injection');
timeAfterDrug = 45;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Firing Rate over Time 5min before and 60 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2); clf;
clear g;
scaling  = sampleRate * secInMin;
downsampleFactor = 250;
timeAfterDrug = 60;
timeBeforeDrug = -5;
x = timeBeforeDrug*scaling+1:downsampleFactor:timeAfterDrug*scaling;
y = downsample(allCellsTimeFRneg5to60minAfterKetamineInjx',downsampleFactor)';

x_min = x./scaling;

g(1,1) = gramm('x',x_min,'y',y); 
g(1,1).stat_summary('type','std');
% g(1,1).stat_smooth();

% g(1,1).geom_line();
g(1,1).axe_property('YLim',[5 10]);
g(1,1).set_title(sprintf('Ketamine-induced FR over Time(%s)',filter), 'FontSize', 40);
g(1,1).set_names('x','Minutes since Injection','y','Firing Rate(Hz)');
g(1,1).set_text_options('base_size',20);
g.draw()
axis square;
set(gcf,'Position',[100 100 1000 1000])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Firing Rate over Time 45 min after injection with KO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf; hold on;
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

plot(smoothdata(nanmean(allCellsTimeFR45minAfterKetamineInjx,1),'gaussian',sampleRate*25),'LineWidth',2,'DisplayName','WT');
plot(smoothdata(nanmean(allCellsTimeFR45minAfterKetamineInjxKO,1),'gaussian',sampleRate*25),'LineWidth',2,'DisplayName','HCN1ko');
% plot(smoothdata(nanmean(allCellsTimeFRcircaKetamineInjx,1),sampleRate*10,'moving'),'LineWidth',2,'DisplayName','Ketamine Injection');
timeAfterDrug = 45;
% modify labels for tick marks
scaling  = sampleRate * secInMin;
tickSteps = 5;
set(gca,'XLim',[0 timeAfterDrug*scaling],'XTick',[0:tickSteps*scaling:timeAfterDrug*scaling])
xticks = get(gca,'xtick');
x = 0:tickSteps:timeAfterDrug;
newlabels = arrayfun(@(x) sprintf('%d', x/scaling), x*scaling, 'un', 0);
set(gca,'xticklabel',newlabels);
title(sprintf('Ketamine-induced FR over Time(WT vs HCN1ko)',filter))
xlabel('Minutes since Injection')
ylabel('Firing Rate (Hz)')
legend;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Firing Rate over Time 45 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot by Mouse

if strcmp(filter, 'WT')
    mice={'G1','G2','G3','G4','G5','HCNd1','HCNe2'};
elseif strcmp(filter, 'KO')
    mice={'HCN1','HCNd2','HCNe1','HCNe3','HCNe4'};
else
    fprintf('Bad filter key. Choose: mec, WT, KO');
end

timeAfterDrug = 45;
scaling  = sampleRate * secInMin;

mouseFRs = nan(numel(mice), size(allCellsTimeFR45minAfterKetamineInjx,2));
for z= 1:numel(mice)
    idx = ~cellfun('isempty',strfind({sessions.name},mice{z}));
    Y = mean(allCellsTimeFR45minAfterKetamineInjx(idx,:));
    interpY = fillmissing(Y,'spline');
    smoothY= sgolayfilt(interpY, 3, scaling+1);
    mouseFRs(z,:) = smoothY;
end
%%
close all; figure(1); clf;
% figure(1); clf; hold on;
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.015 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
% hold on;


frPlot.x = 0:1/scaling:timeAfterDrug;
frPlot.y = mouseFRs;
% Y = allCellsTimeFR45minAfterKetamineInjx;
% interpY = fillmissing(Y,'spline');

% frPlot.y = interpY;





clear g;

x = frPlot.x;
y = frPlot.y;

g(1,1) = gramm('x',x,'y',y); 
g(1,1).stat_summary('type','std');
% g(1,1).stat_smooth();
% g(1,1).geom_line();
g(1,1).axe_property('YLim',[0 15]);
g(1,1).set_title("Ketamine's Effect on Firing Rate by Time", 'FontSize', 30);
g(1,1).set_names('x','Minutes','y','Firing Rate(Hz)');
g(1,1).set_text_options('base_size',30);

g.draw()


%%
fs= 50;
pwelch(allCellsTimeFR60minAfterKetamineInjx(5,:), 10*fs, [] , [], fs)
pwelch(allCellsTimeFRcircaKetamineInjx(5,:), 10*fs, [] , [], fs)