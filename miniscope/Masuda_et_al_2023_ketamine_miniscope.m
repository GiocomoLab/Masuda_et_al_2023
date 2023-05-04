%% This code is used to generate the figure related to miniscope imaging in Masuda et al., 2023
%% Calculate and store information score in batch by comboSession
load ('neuronIndividualsf.mat');
load ('behavIndividualsf.mat');
load ('thresh.mat');
meanFRindividuals = cell(1,length(neuronIndividualsf));
infoPerSecIndividuals = cell(1,length(neuronIndividualsf));
infoPerSpikeIndividuals = cell(1,length(neuronIndividualsf));
firingrateAll = cell(1,length(neuronIndividualsf));
countAll = cell(1,length(neuronIndividualsf));
countTime = cell(1,length(neuronIndividualsf));

% Calculate mean firing rate for each individual sessions
for k = 1:length(neuronIndividualsf)
    [firingrateAll{k},countAll{k},countTime{k}] = calculate_firing_ratemap(neuronIndividualsf{k},behavIndividualsf{k},thresh,3);
    for m = 1:size(firingrateAll{1,k},2)
       meanFRindividuals{1,k}(m,1)= sum(sum(countAll{1,k}{1,m}))/sum(sum(countTime{1,k}));
    end
end
save firingrateAll.mat firingrateAll countAll countTime meanFRindividuals;

%% Incorporate the firing rate information and then plot the final spatial coding heat map
plot_rate_map_longitudinal(neuronIndividualsf,behavIndividualsf,firingrateAll,placecell.all(101:150),thresh,'S',2); %1:size(neuronIndividualsf{1}.C,1); PlaceCells_batch{3}(100:end,1)
% Print figure into vector eps files
filepath = pwd;
neuronSelected = 169;
fig = openfig(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected), '_ratemap.fig']));
print (fig, '-painters', '-depsc', fullfile(filepath,'RatemapFigures',['Cell',num2str(neuronSelected),'_ratemap.eps']));
% open a group of plots;
neuronSelected = dispcm.p(1:30);
for ii = 1:length(neuronSelected)
    if exist(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected(ii)), '_ratemap.fig']), 'file')
        openfig(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected(ii)), '_ratemap.fig']));
    end
end
%% Identify place cells by shuffling spikes for each mouse
datapath = {};
for n = 1:length(datapath)
    cd(datapath{n});
    load('sessions.mat');
    load('neuronIndividualsf.mat');
    load('behavIndividualsf.mat');
    load('thresh.mat');
    PlaceCells_batch = cell(1,length(sessions));
    parfor ii = 1:length(sessions)
        filepath = [pwd '\' sessions{ii}]; %ii
        neuron0 = neuronIndividualsf{ii}; %ii
        behav = behavIndividualsf{ii}; %ii
        [PlaceCells_both,Tinfo,SIthresh] = permutingSpike_speed(neuron0,behav,thresh);
        parsave(filepath, 'PlaceCell_Tinfo.mat', PlaceCells_both, Tinfo);
        PlaceCells_batch{ii} = PlaceCells_both{1}; %defined by bit/sec
    end
    save PlaceCells_batch.mat PlaceCells_batch;
    % redefine place cells by filtering out low fr neurons
    load ('firingrateAll.mat', 'meanFRindividuals');
    % get defined place cells in baseline session
    placecell = struct;
    placecell.highfrCell = find(meanFRindividuals{2} > 0.1); 
    placecell.all = intersect(placecell.highfrCell, PlaceCells_batch{2});
    save placecell.mat placecell;
    % get place cells for all the sessions by filtering out low fr neurons
    highfrcell = cellfun(@(x) find(x > 0.1), meanFRindividuals, 'UniformOutput', false); % filter out low fr neurons
    placecells_batch = cellfun(@intersect, PlaceCells_batch, highfrcell, 'Uniformoutput', false);
    save('PlaceCells_batch.mat', 'placecells_batch', '-append');
end
%% Compute rate map correlation for each mouse
datapath = {};
for n = 1:length(datapath)
    cd(datapath{n});
    load('neuronIndividualsf.mat','neuronIndividualsf');
    load('behavIndividualsf.mat','behavIndividualsf');
    load('thresh.mat');
    load('placecell.mat')
    niter = 20;
    corrX = cell(length(neuronIndividualsf),length(neuronIndividualsf),niter); 
    corrX_avg = NaN(length(neuronIndividualsf),length(neuronIndividualsf),niter);
    for ii = 1:niter
        [ratemapt_Xm, frXm, countXm, timeXm] = match_smooth_trim_ratemap(...
            neuronIndividualsf, behavIndividualsf, thresh, 'xsession');
        [corr_XLm, avg_corr_XLm] = get_corr_matrix([], ratemapt_Xm, [], 'match');
        corrX(:,:,ii) = corr_XLm;
        corrX_avg(:,:,ii) = avg_corr_XLm;
    end
    corrX_pc = cellfun(@(x) x(placecell.all,:), corrX, 'uni', false);
    corrX_pcavg = cellfun(@mean, corrX_pc);
    save('Corr_match.mat', 'corrX', 'corrX_avg', 'corrX_pc', 'corrX_pcavg','-v7.3');
end
% average iterated correlations
datapath = {};
for n = 1:length(datapath)
    cd(datapath{n});
    corr_pc = {};
    load('Corr_match.mat')
    for ii = 1:size(corrX_pc,3)
    corr12 = cell2mat(corrX_pc(1,2,:));
    corr13 = cell2mat(corrX_pc(1,3,:));
    corr14 = cell2mat(corrX_pc(1,4,:));
    corr23 = cell2mat(corrX_pc(2,3,:));
    corr24 = cell2mat(corrX_pc(2,4,:));
    corr34 = cell2mat(corrX_pc(3,4,:));
    end
    corr_pc{1,2} = mean(corr12,3);
    corr_pc{1,3} = mean(corr13,3);
    corr_pc{1,4} = mean(corr14,3);
    corr_pc{2,3} = mean(corr23,3);
    corr_pc{2,4} = mean(corr24,3);
    corr_pc{3,4} = mean(corr34,3);
    save('Corr_match.mat', 'corr_pc', '-append')
end

%% plot correlation matrix for all mice (panel i)
load('E:\Miniscope_Ketamine_analysis_080519\Miniscope_Ketamine_results\Data_location.mat')
figure
for n = 1:length(datapath)
    cd(datapath{n});
    load('Corr_match.mat')
    subplot(2,4,n);
    imagesc(mean(corrX_pcavg,3))
    axis image
    title(['mouse ',num2str(n)])
end

%% plot correlation histogram (panel g-h)
load('E:\Miniscope_Ketamine_analysis_080519\Miniscope_Ketamine_results\Data_location.mat')
corr12 = [];
corr23 = [];
corr24 = [];
for n = 1:length(datapath)
    cd(datapath{n});
    load('Corr_match.mat', 'corr_pc')
    corr12 = [corr12; corr_pc{1,2}];
    corr23 = [corr23; corr_pc{2,3}];
    corr24 = [corr24; corr_pc{2,4}];
end
corr12_avg = [];
corr23_avg = [];
corr24_avg = [];
for n = 1:length(datapath)
    cd(datapath{n});
    load('Corr_match.mat', 'corr_pc')
    corr12_avg = [corr12_avg; mean(corr_pc{1,2})];
    corr23_avg = [corr23_avg; mean(corr_pc{2,3})];
    corr24_avg = [corr24_avg; mean(corr_pc{2,4})];
end

% plot panel g
figure;
histogram(corr12)
hold on
histogram(corr23)
% histogram(corr24)
xlabel('Correlation')
ylabel('Num of place cells')
legend({'Bsl1 vs. Bsl2', 'Bsl2 vs. Ketamine'})
[p,~] = signrank(corr12_avg, corr23_avg);
figure;
plot([corr12_avg,corr23_avg]', 'o-', 'color','k')
xlim([0.5,2.5])
xticks([1:2]);
xticklabels({'bsl2 vs. bsl1', 'bsl2 vs. K'})
ylabel('Spatial correlation')
text(1.5,0.65, ['P=', num2str(p)])

%plot panel h
figure;
histogram(corr12)
hold on
histogram(corr24, 'FaceColor', 'g', 'FaceAlpha', 0.3)
xlabel('Correlation')
ylabel('Num of place cells')
legend({'Bsl1 vs. Bsl2', 'Bsl2 vs. Post-Bsl'})
[p,~] = signrank(corr12_avg, corr24_avg);
figure;
plot([corr12_avg,corr24_avg]', 'o-', 'color','k')
xlim([0.5,2.5])
xticks([1:2]);
xticklabels({'bsl2 vs. bsl1', 'bsl2 vs. p-bsl'})
ylabel('Spatial correlation')
text(1.5,0.68, ['P=', num2str(p)])

%% num of place cells (panel f)
load('E:\Miniscope_Ketamine_analysis_080519\Miniscope_Ketamine_results\Data_location.mat')
numpc = [];
for n = 1:length(datapath)
    cd(datapath{n});
    load('PlaceCells_batch.mat', 'placecells_batch')
    x = cellfun(@numel, placecells_batch, 'uni', false);
    numpc = [numpc; x];
end
% plot panel f
numpc = cell2mat(numpc);
[~,p1] = ttest(numpc(:,2),numpc(:,1));
[~,p2] = ttest(numpc(:,2),numpc(:,3));
[~,p3] = ttest(numpc(:,2),numpc(:,4));
figure;
plot(numpc', '-o', 'Color', 'k');
xlim([0.5 4.5]);
xticks([1 2 3 4])
xticklabels({'bsl1','bsl2','K','p-bsl'})
ylabel('Number of place cells')
text(1,450, ['P=',num2str(p1)])
text(3,450, ['P=',num2str(p2)])
text(4,450, ['P=',num2str(p3)])

%% calculate metrics by using pair-wise matched data (for panel c-e)
load('E:\Miniscope_Ketamine_analysis_080519\Miniscope_Ketamine_results\Data_location.mat')
for n = 1:length(datapath)
    cd(datapath{n});
    load('neuronIndividualsf.mat','neuronIndividualsf');
    load('behavIndividualsf.mat','behavIndividualsf');
    load('thresh.mat');
    niter = 20;
    bpsecXm = NaN(length(thresh),6,niter); 
    bpspikeXm = NaN(length(thresh),6,niter); 
    meanfrXm = NaN(length(thresh),6,niter); 
    peakfrXm = NaN(length(thresh),6,niter);
for iter = 1:niter
        [ratemapt_Xm, frXm, countXm, timeXm] = match_smooth_trim_ratemap(...
            neuronIndividualsf, behavIndividualsf, thresh, 'xsession');       
        %b1-b2
        meanfr = cell2mat(cellfun(@(x) sum(sum(x)),countXm{1,2},'uni',0))'/sum(sum(timeXm{1,2}));
        for jj = 1:length(meanfr)
            [bpsecXm(jj,1,iter), bpspikeXm(jj,1,iter)] = Doug_spatialInfo(...
                frXm{1,2}{jj}, meanfr(jj), timeXm{1,2}, 0.1);
        end
        meanfrXm(:,1,iter) = meanfr;
        peakfrXm(:,1,iter) = cell2mat(cellfun(@(x) max(max(x)), ratemapt_Xm{1,2},'uni',0))';
        
        meanfr = cell2mat(cellfun(@(x) sum(sum(x)),countXm{2,1},'uni',0))'/sum(sum(timeXm{2,1}));
        for jj = 1:length(meanfr)
            [bpsecXm(jj,2,iter), bpspikeXm(jj,2,iter)] = Doug_spatialInfo(...
                frXm{2,1}{jj}, meanfr(jj), timeXm{2,1}, 0.1);
        end
        meanfrXm(:,2,iter) = meanfr;
        peakfrXm(:,2,iter) = cell2mat(cellfun(@(x) max(max(x)), ratemapt_Xm{2,1},'uni',0))';

        %b2-k/b2-b3
        meanfr = cell2mat(cellfun(@(x) sum(sum(x)),countXm{2,3},'uni',0))'/sum(sum(timeXm{2,3}));
        for jj = 1:length(meanfr)
        [bpsecXm(jj,3,iter), bpspikeXm(jj,3,iter)] = Doug_spatialInfo(...
            frXm{2,3}{jj}, meanfr(jj), timeXm{2,3}, 0.1);
        end
        meanfrXm(:,3,iter) = meanfr;
        peakfrXm(:,3,iter) = cell2mat(cellfun(@(x) max(max(x)), ratemapt_Xm{2,3},'uni',0))';
        
        meanfr = cell2mat(cellfun(@(x) sum(sum(x)),countXm{3,2},'uni',0))'/sum(sum(timeXm{3,2}));
        for jj = 1:length(meanfr)
            [bpsecXm(jj,4,iter), bpspikeXm(jj,4,iter)] = Doug_spatialInfo(...
                frXm{3,2}{jj}, meanfr(jj), timeXm{3,2}, 0.1);
        end
        meanfrXm(:,4,iter) = meanfr;
        peakfrXm(:,4,iter) = cell2mat(cellfun(@(x) max(max(x)), ratemapt_Xm{3,2},'uni',0))';

        %b2-pb/b2-b4
        meanfr = cell2mat(cellfun(@(x) sum(sum(x)),countXm{2,4},'uni',0))'/sum(sum(timeXm{2,4}));
        for jj = 1:length(meanfr)
            [bpsecXm(jj,5,iter), bpspikeXm(jj,5,iter)] = Doug_spatialInfo(...
                frXm{2,4}{jj}, meanfr(jj), timeXm{2,4}, 0.1);
        end
        meanfrXm(:,5,iter) = meanfr;
        peakfrXm(:,5,iter) = cell2mat(cellfun(@(x) max(max(x)), ratemapt_Xm{2,4},'uni',0))';
        
        meanfr = cell2mat(cellfun(@(x) sum(sum(x)),countXm{4,2},'uni',0))'/sum(sum(timeXm{4,2}));
        for jj = 1:length(meanfr)
            [bpsecXm(jj,6,iter), bpspikeXm(jj,6,iter)] = Doug_spatialInfo(...
                frXm{4,2}{jj}, meanfr(jj), timeXm{4,2}, 0.1);
        end
        meanfrXm(:,6,iter) = meanfr;
        peakfrXm(:,6,iter) = cell2mat(cellfun(@(x) max(max(x)), ratemapt_Xm{4,2},'uni',0))';
end
    %final results
    metricsXm = struct;
    metricsXm.meanfr = mean(meanfrXm,3);
    metricsXm.peakfr = mean(peakfrXm,3);
    metricsXm.bpsec = mean(bpsecXm,3);
    metricsXm.bpspike = mean(bpspikeXm,3);
    save('Metrics_match.mat', 'metricsXm');
end

%% For plotting figure panel c-e
metrics = struct;
metrics.meanfr = [];
metrics.peakfr = [];
metrics.bpspike = [];
metrics.bpsec = [];
metrics.meanfr_avg = [];
metrics.peakfr_avg = [];
metrics.bpspike_avg = [];
metrics.bpsec_avg = [];
metrics.mouse = [];
for n = 1:length(datapath)
    cd(datapath{n});
    load('Metrics_match.mat')
    metrics.meanfr = [metrics.meanfr; metricsXm.meanfr];
    metrics.peakfr = [metrics.peakfr; metricsXm.peakfr];
    metrics.bpspike = [metrics.bpspike; metricsXm.bpspike];
    metrics.bpsec = [metrics.bpsec; metricsXm.bpsec];
    metrics.meanfr_avg = [metrics.meanfr_avg; median(metricsXm.meanfr)];
    metrics.peakfr_avg = [metrics.peakfr_avg; median(metricsXm.peakfr)];
    metrics.bpspike_avg = [metrics.bpspike_avg; median(metricsXm.bpspike)];
    metrics.bpsec_avg = [metrics.bpsec_avg; median(metricsXm.bpsec)];
    metrics.mouse = [metrics.mouse; repelem(n, size(metricsXm.meanfr,1))'];
end
save('E:\Miniscope_Ketamine_analysis_080519\Miniscope_Ketamine_results\metricsXm.mat', 'metrics')

%% plot panel c
data1 = metrics.meanfr(:,1:2);
data_avg1 = metrics.meanfr_avg(:,1:2);
[p1,~] = signrank(data_avg1(:,1),data_avg1(:,2));
r_low1 = abs(quantile(data1,0.25)- nanmedian(data1));
r_high1 = abs(quantile(data1,0.75)-nanmedian(data1));
data2 = metrics.meanfr(:,3:4);
data_avg2 = metrics.meanfr_avg(:,3:4);
[p2,~] = signrank(data_avg2(:,1),data_avg2(:,2));
r_low2 = abs(quantile(data2,0.25)- nanmedian(data2));
r_high2 = abs(quantile(data2,0.75)-nanmedian(data2));
data3 = metrics.meanfr(:,5:6);
data_avg3 = metrics.meanfr_avg(:,5:6);
[p3,~] = signrank(data_avg3(:,1),data_avg3(:,2));
r_low3 = abs(quantile(data3,0.25)- nanmedian(data3));
r_high3 = abs(quantile(data3,0.75)-nanmedian(data3));
yl = 'mean Ca2+ event rate';
ylm = [0.14,0.36];

ax = figure;
set(ax, 'Position', [200, 200, 800, 600]);
subplot (1,3,1)
errorbar([0.9:1.9], nanmedian(data1),...
    r_low1, r_high1, '.-', 'MarkerSize', 30, 'Color', 'k');
hold on
plot(data_avg1',':o','color',[109 110 113]/255)
xlim([0.3,2.5])
xticks([1:2]);
xticklabels({'bsl1', 'bsl2'});
ylabel(yl)
text(1.5, 3.2, ['P1=',num2str(p1)])
axis square
ylim(ylm)

subplot (1,3,2)
errorbar([0.9:1.9], nanmedian(data2),...
    r_low2, r_high2, '.-', 'MarkerSize', 30, 'Color', 'k');
hold on
plot(data_avg2',':o','color',[109 110 113]/255)
xlim([0.3,2.5])
xticks([1:2]);
xticklabels({'bsl2', 'K'});
ylabel(yl)
text(1.5, 3.2, ['P2=',num2str(p2)])
axis square
ylim(ylm)

subplot (1,3,3)
errorbar([0.9:1.9], nanmedian(data3),...
    r_low3, r_high3, '.-', 'MarkerSize', 30, 'Color', 'k');
hold on
plot(data_avg3',':o','color',[109 110 113]/255)
xlim([0.3,2.5])
xticks([1:2]);
xticklabels({'bsl2', 'p-bsl'});
ylabel(yl)
text(1.5, 3.2, ['P3=',num2str(p3)])
axis square
ylim(ylm)

%% plot panel d
data1 = metrics.peakfr(:,1:2);
data_avg1 = metrics.peakfr_avg(:,1:2);
[p1,~] = signrank(data_avg1(:,1),data_avg1(:,2));
r_low1 = abs(quantile(data1,0.25)- nanmedian(data1));
r_high1 = abs(quantile(data1,0.75)-nanmedian(data1));
data2 = metrics.peakfr(:,3:4);
data_avg2 = metrics.peakfr_avg(:,3:4);
[p2,~] = signrank(data_avg2(:,1),data_avg2(:,2));
r_low2 = abs(quantile(data2,0.25)- nanmedian(data2));
r_high2 = abs(quantile(data2,0.75)-nanmedian(data2));
data3 = metrics.peakfr(:,5:6);
data_avg3 = metrics.peakfr_avg(:,5:6);
[p3,~] = signrank(data_avg3(:,1),data_avg3(:,2));
r_low3 = abs(quantile(data3,0.25)- nanmedian(data3));
r_high3 = abs(quantile(data3,0.75)-nanmedian(data3));
yl = 'peak Ca2+ event rate';
ylm = [0.6,2];

ax = figure;
set(ax, 'Position', [200, 200, 800, 600]);
subplot (1,3,1)
errorbar([0.9:1.9], nanmedian(data1),...
    r_low1, r_high1, '.-', 'MarkerSize', 30, 'Color', 'k');
hold on
plot(data_avg1',':o','color',[109 110 113]/255)
xlim([0.3,2.5])
xticks([1:2]);
xticklabels({'bsl1', 'bsl2'});
ylabel(yl)
text(1.5, 3.2, ['P1=',num2str(p1)])
axis square
ylim(ylm)

subplot (1,3,2)
errorbar([0.9:1.9], nanmedian(data2),...
    r_low2, r_high2, '.-', 'MarkerSize', 30, 'Color', 'k');
hold on
plot(data_avg2',':o','color',[109 110 113]/255)
xlim([0.3,2.5])
xticks([1:2]);
xticklabels({'bsl2', 'K'});
ylabel(yl)
text(1.5, 3.2, ['P2=',num2str(p2)])
axis square
ylim(ylm)

subplot (1,3,3)
errorbar([0.9:1.9], nanmedian(data3),...
    r_low3, r_high3, '.-', 'MarkerSize', 30, 'Color', 'k');
hold on
plot(data_avg3',':o','color',[109 110 113]/255)
xlim([0.3,2.5])
xticks([1:2]);
xticklabels({'bsl2', 'p-bsl'});
ylabel(yl)
text(1.5, 3.2, ['P3=',num2str(p3)])
axis square
ylim(ylm)

%% plot panel e
data1 = metrics.bpsec(:,1:2);
data_avg1 = metrics.bpsec_avg(:,1:2);
[p1,~] = signrank(data_avg1(:,1),data_avg1(:,2));
r_low1 = abs(quantile(data1,0.25)- nanmedian(data1));
r_high1 = abs(quantile(data1,0.75)-nanmedian(data1));
data2 = metrics.bpsec(:,3:4);
data_avg2 = metrics.bpsec_avg(:,3:4);
[p2,~] = signrank(data_avg2(:,1),data_avg2(:,2));
r_low2 = abs(quantile(data2,0.25)- nanmedian(data2));
r_high2 = abs(quantile(data2,0.75)-nanmedian(data2));
data3 = metrics.bpsec(:,5:6);
data_avg3 = metrics.bpsec_avg(:,5:6);
[p3,~] = signrank(data_avg3(:,1),data_avg3(:,2));
r_low3 = abs(quantile(data3,0.25)- nanmedian(data3));
r_high3 = abs(quantile(data3,0.75)-nanmedian(data3));
yl = 'bits / sec';
ylm = [0.4,0.8];

ax = figure;
set(ax, 'Position', [200, 200, 800, 600]);
subplot (1,3,1)
errorbar([0.9:1.9], nanmedian(data1),...
    r_low1, r_high1, '.-', 'MarkerSize', 30, 'Color', 'k');
hold on
plot(data_avg1',':o','color',[109 110 113]/255)
xlim([0.3,2.5])
xticks([1:2]);
xticklabels({'bsl1', 'bsl2'});
ylabel(yl)
text(1.5, 3.2, ['P1=',num2str(p1)])
axis square
ylim(ylm)

subplot (1,3,2)
errorbar([0.9:1.9], nanmedian(data2),...
    r_low2, r_high2, '.-', 'MarkerSize', 30, 'Color', 'k');
hold on
plot(data_avg2',':o','color',[109 110 113]/255)
xlim([0.3,2.5])
xticks([1:2]);
xticklabels({'bsl2', 'K'});
ylabel(yl)
text(1.5, 3.2, ['P2=',num2str(p2)])
axis square
ylim(ylm)

subplot (1,3,3)
errorbar([0.9:1.9], nanmedian(data3),...
    r_low3, r_high3, '.-', 'MarkerSize', 30, 'Color', 'k');
hold on
plot(data_avg3',':o','color',[109 110 113]/255)
xlim([0.3,2.5])
xticks([1:2]);
xticklabels({'bsl2', 'p-bsl'});
ylabel(yl)
text(1.5, 3.2, ['P3=',num2str(p3)])
axis square
ylim(ylm)

