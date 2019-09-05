%% Load data, get spikes from good cells

addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));

load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/G4_190620_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');
load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/HCN1_190619_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');
load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/G3_190708_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');
cells_to_plot = sp.cids(sp.cgs==2); 
nCells = size(cells_to_plot, 2);

% compute some useful information (like spike depths)
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

% get spike depths
spike_depth = nan(numel(cells_to_plot),1);
for k = 1:numel(cells_to_plot)
    spike_depth(k) = median(spikeDepths(sp.clu==cells_to_plot(k)));
end

% sort cells_to_plot by spike_depth (descending)
[spike_depth,sort_idx] = sort(spike_depth,'descend');
cells_to_plot = cells_to_plot(sort_idx);
%%
% 
% % get spike times and index into post for cell k 
% spike_t = sp.st(sp.clu==cells_to_plot(k));
% 
% [~,~,spike_idx] = histcounts(spike_t,post);
% [~,~,spikeTrial_idx] = histcounts(spike_t,trial); % NOT SURE IF WE STILL NEED THIS
% 
% all_fr = [];
% ds_factor = 1;
% smoothSigma = 10;
% 
% % for cell k, iteratively calculate the firing rate for each trial
% for i = 1:max(trial)
%     fr = calcSmoothedFR_Time(post, spike_t, ds_factor, smoothSigma); 
%     itrial_kfr = calcSmoothedFR_SpatialBin(spike_idx(i==trial(spike_idx)), posx, p, trackEnd);
%     singleCellallTrialsFR(i,:) = itrial_kfr;
% %         figure(1)
% %         plot(itrial_kfr); %plot FR for single trial in single cell
% %         figure(2) % plot rasterplot for single trial in single cell
% %         plot(posx(spike_idx(i==trial(spike_idx))),trial(spike_idx(i==trial(spike_idx))),'k.');      
% end



%% Calculating firing rate

% calculate firing rate by time
end_recording = max(post);
all_fr = [];
ds_factor = 50;
smoothSigma = 10;
trial_ds = downsample(trial, ds_factor); 
for i = 1:nCells
   
    % get spike times for cell i
    spike_t = sp.st(sp.clu==cells_to_plot(i));
    
    % calculate firing rate
    fr = calcSmoothedFR_Time(post, spike_t, ds_factor, smoothSigma);
%     fr = calcSmoothedFR_SpatialBin(idx, posx, p, TrackEnd);
    all_fr(i, :) = fr;
        
end

all_fr = all_fr';

figure(1);
plot(mean(all_fr,2));

% close all;
% clear g;
% h = figure('Position',[100 100 1800 600]);
% g(1,1) = gramm('x',trial,'y',mean(all_fr,2));
% g(1,1).geom_line;
% g(1,1).set_names('x','Lick Position in VR','y','trial number');
% g.draw

%% PCA 



%[eigvecs_all, transformed_all, eigvals_all, ~, ~, mu_all] = pca(all_fr);

[coeff, score, latent, tsquared, explained] = pca(all_fr);
figure(1);
plot(explained(1:10));
title('Variance Explained')
% Screen plot to visually determine how many principal components are
% needed to capture variance in the data
for i = 1:8
   sum(explained(1:i)) 
end

%%

% use the elbow in the plot to determinee how many principal components are
% needed. 

nComponents = 3;
% Project onto the first three PCA components and plot as an animated
% line

firstThreeComponentsAll = coeff(:, [1:3]);

projectedDataAll = all_fr * firstThreeComponentsAll;

x_all=projectedDataAll(:,1);
y_all=projectedDataAll(:,2);
z_all=projectedDataAll(:,3);
%%
close all;
figure(2); 
plot3(x_all(trial_ds<51),y_all(trial_ds<51),z_all(trial_ds<51),'b.'); hold on;
plot3(x_all(trial_ds>50 & trial_ds<101),y_all(trial_ds>50 & trial_ds<101),z_all(trial_ds>50 & trial_ds<101),'k.');
plot3(x_all(trial_ds>100 & trial_ds<151),y_all(trial_ds>100 & trial_ds<151),z_all(trial_ds>100& trial_ds<151),'r.');
plot3(x_all(trial_ds>150),y_all(trial_ds>150),z_all(trial_ds>150),'g.'); hold on;
legend('Baseline','ControlInjx','Ketamine:1-50','Ketamine:51-200');


figure(3); 
plot(x_all(trial_ds<51),y_all(trial_ds<51),'color','b'); hold on;
plot(x_all(trial_ds>50 & trial_ds<101),y_all(trial_ds>50 & trial_ds<101),'color','k');
plot(x_all(trial_ds>100 & trial_ds<151),y_all(trial_ds>100 & trial_ds<151),'color','r');
plot(x_all(trial_ds>150),y_all(trial_ds>150),'color','g'); hold on;
legend('Baseline','ControlInjx','Ketamine:1-50','Ketamine:51-200');
%% Variance explained by PCA on pre-, early, mid, late, and post-VR
% X*v*v' 

Data_Indices = {[1:pre_VR_frames], [pre_VR_frames+1: 2* pre_VR_frames], ...
    [ceil(size(all_fr, 1)/2) + 1: ceil(size(all_fr, 1)/2) + pre_VR_frames], ...
    [size(all_fr, 1)-2*pre_VR_frames+1:size(all_fr, 1)-pre_VR_frames], ...
    [size(all_fr, 1)-pre_VR_frames+1:size(all_fr, 1)], ...
    [pre_VR_frames+1:size(all_fr, 1)-pre_VR_frames]};

labels = {('Pre-VR'), ('Early-VR'), ('Mid-VR'), ('Late-VR'), ('Post-VR'), ...
    ('All of VR')};

titles = {['Variance Explained by Pre-VR PCs'], ['Variance Explained by Early-VR PCs'], ...
    ['Variance Explained by Mid-VR PCs'], ['Variance Explained by Late-VR PCs'], ...
    ['Variance Explained by Post-VR PCs'], ['Variance Explained by All VR PCs']};

% Populate PCA struct
PCA_st.Data_Indices = Data_Indices;
PCA_st.labels = labels;
PCA_st.titles = titles;
PCA_st.var_expl = {};
PCA_st.projected = {};

for i = 1:size(PCA_st.Data_Indices, 2)
    
    dat = all_fr(PCA_st.Data_Indices{i}, :); % data on which to do PCA
    [eigvecs, transformed, eigvals] = pca(dat); 
    
    features = eigvecs(:, [1:nComponents]); % get eigenvectors from top n eigenvalues
    PCA_st.projected{i} = VR_fr * features;
    
    xhat = VR_fr * features * features'; 
    PCA_st.var_expl{i} = 1 - var(VR_fr-xhat)/var(VR_fr);
    
end

bar(cell2mat(PCA_st.var_expl));
set(gca, 'XTickLabel', PCA_st.labels, 'FontSize', 10);
text(1:length(PCA_st.var_expl), cell2mat(PCA_st.var_expl), num2str(cell2mat(PCA_st.var_expl')),'vert','bottom','horiz','center'); 
box off
ylabel('Variance Explained', 'FontSize', 20);
title('Variance Explained vs. PCs from Different Portions of Data', 'FontSize', 18)

%% PLOTS 

% PCA ALL
figure()
plot3(x_all(1:pre_VR_frames), y_all(1:pre_VR_frames), z_all(1:pre_VR_frames), 'r');
hold on
plot3(x_all(pre_VR_frames+1:length(x_all)-pre_VR_frames-pre_VR_frames), ...
    y_all(pre_VR_frames+1:length(x_all)-pre_VR_frames-pre_VR_frames), ...
    z_all(pre_VR_frames+1:length(x_all)-pre_VR_frames-pre_VR_frames), 'b');
hold on
plot3(x_all(length(x_all)-pre_VR_frames-pre_VR_frames:length(x_all)), ...
    y_all(length(x_all)-pre_VR_frames-pre_VR_frames:length(x_all)), ...
    z_all(length(x_all)-pre_VR_frames-pre_VR_frames:length(x_all)), 'm');
legend('preVR', 'VR', 'postVR');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('PCA_{all} before, during, and after VR');

figure();
plot(x_all(1:pre_VR_frames), y_all(1:pre_VR_frames), 'r');
hold on
plot(x_all(pre_VR_frames+1:length(x_all)-pre_VR_frames-pre_VR_frames), ...
    y_all(pre_VR_frames+1:length(x_all)-pre_VR_frames-pre_VR_frames), 'b');
hold on
plot(x_all(length(x_all)-pre_VR_frames-pre_VR_frames:length(x_all)), ...
    y_all(length(x_all)-pre_VR_frames-pre_VR_frames:length(x_all)), 'm');
legend('preVR', 'VR', 'postVR');
xlabel('PC1');
ylabel('PC2');
title('PCA_{all} before, during, and after VR');


% PCA_VR
figure();
plot3(x_vr(1:pre_VR_frames), y_vr(1:pre_VR_frames), z_vr(1:pre_VR_frames), 'r');
hold on
plot3(x_vr(pre_VR_frames+1:length(x_vr)-pre_VR_frames-pre_VR_frames), ...
    y_vr(pre_VR_frames+1:length(x_vr)-pre_VR_frames-pre_VR_frames), ...
    z_vr(pre_VR_frames+1:length(x_vr)-pre_VR_frames-pre_VR_frames), 'b');
hold on
plot3(x_vr(length(x_vr)-pre_VR_frames-pre_VR_frames:length(x_vr)), ...
    y_vr(length(x_vr)-pre_VR_frames-pre_VR_frames:length(x_vr)), ...
    z_vr(length(x_vr)-pre_VR_frames-pre_VR_frames:length(x_vr)), 'm');
legend('preVR', 'VR', 'postVR');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('PCA_{VR} before, during, and after VR');

figure();
plot(x_vr(1:pre_VR_frames), y_vr(1:pre_VR_frames), 'r');
hold on
plot(x_vr(pre_VR_frames+1:length(x_vr)-pre_VR_frames-pre_VR_frames), ...
    y_vr(pre_VR_frames+1:length(x_vr)-pre_VR_frames-pre_VR_frames), 'b');
hold on
plot(x_vr(length(x_vr)-pre_VR_frames-pre_VR_frames:length(x_vr)), ...
    y_vr(length(x_vr)-pre_VR_frames-pre_VR_frames:length(x_vr)), 'm');
legend('preVR', 'VR', 'postVR');
xlabel('PC1');
ylabel('PC2');
title('PCA_{VR} before, during, and after VR');


% PCA_PRE
figure()
plot3(x_pre(1:pre_VR_frames), y_pre(1:pre_VR_frames), z_pre(1:pre_VR_frames), 'r');
hold on
plot3(x_pre(pre_VR_frames+1:length(x_pre)-pre_VR_frames-pre_VR_frames), ...
    y_pre(pre_VR_frames+1:length(x_pre)-pre_VR_frames-pre_VR_frames), ...
    z_pre(pre_VR_frames+1:length(x_pre)-pre_VR_frames-pre_VR_frames), 'b');
hold on
plot3(x_pre(length(x_pre)-pre_VR_frames-pre_VR_frames:length(x_pre)), ...
    y_pre(length(x_pre)-pre_VR_frames-pre_VR_frames:length(x_pre)), ...
    z_pre(length(x_pre)-pre_VR_frames-pre_VR_frames:length(x_pre)), 'm');
legend('preVR', 'VR', 'postVR');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('PCA_{pre} before, during, and after VR');

figure();
plot(x_pre(1:pre_VR_frames), y_pre(1:pre_VR_frames), 'r');
hold on
plot(x_pre(pre_VR_frames+1:length(x_pre)-pre_VR_frames-pre_VR_frames), ...
    y_pre(pre_VR_frames+1:length(x_pre)-pre_VR_frames-pre_VR_frames), 'b');
hold on
plot(x_pre(length(x_pre)-pre_VR_frames-pre_VR_frames:length(x_pre)), ...
    y_pre(length(x_pre)-pre_VR_frames-pre_VR_frames:length(x_pre)), 'm');
legend('preVR', 'VR', 'postVR');
xlabel('PC1');
ylabel('PC2');
title('PCA_{pre} before, during, and after VR');




