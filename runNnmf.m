function doRICA(matPath)

addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));
% 
% load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/G4_190620_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');
% load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/HCN1_190619_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');
% load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/G3_190708_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');

load(matPath);

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

%% Calculating firing rate

% calculate firing rate by time
all_fr = []; s 
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

%% Z score 
imgDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/img';
[~,filename] = fileparts(matPath);
%plot(mean(all_fr)
zFR = zscore(all_fr);

%% NNMF
figure(1);
k = 5;
clf;
[w,h] = nnmf(zFR,k);
figure(1);
for i = 1:k
    plot(trial_ds,w(:,i)); hold on;
end

figure(2); hold on;
clf;
sortH = sortrows(h',1);
imagesc(sortH);
colorbar;

figure(3)
hLines = gscatter(w(:,1),w(:,2));
hold on; biplot(max(w(:))*h','VarLabels',{'sl' 'sw' 'pl' 'pw'},'positive',true); hold off;
axis([0 12 0 12]);
legend(hLines)

%% PCA

[coeff, score, latent, ~, explained] = pca(zFR);
clf;
figure(1);
set(0,'defaultAxesFontSize',25)
set(gcf,'Position',[100 100 2000 1000]); 
subplot(2, 2, 2);
plot(cumsum(explained(1:10)));
title('% Variance Explained');
xlabel('Principal Components');
ylabel('% Variance Explained');

% Project onto the first three PCA components and plot as an animated
% line
firstThreeComponentsAll = coeff(:, [1:3]);

projectedDataAll = all_fr * firstThreeComponentsAll;

x_all=projectedDataAll(:,1);
y_all=projectedDataAll(:,2);
z_all=projectedDataAll(:,3);

trialBlocks = nan(size(trial_ds,1),3);
trialBlocks(trial_ds <= 50,1) = 0; trialBlocks(trial_ds <= 50,2) = 0; trialBlocks(trial_ds <= 50,3) = 1;
trialBlocks(trial_ds>=51 & trial_ds<=100,1) = 0; trialBlocks(trial_ds>=51 & trial_ds<=100,2) = 0; trialBlocks(trial_ds>=51 & trial_ds<=100,3) = 0;
trialBlocks(trial_ds>=101 & trial_ds<=150,1) = 1; trialBlocks(trial_ds>=101 & trial_ds<=150,2) = 0; trialBlocks(trial_ds>=101 & trial_ds<=150,3) = 0;
trialBlocks(trial_ds >= 151,1) = 0; trialBlocks(trial_ds >= 151,2) = 1; trialBlocks(trial_ds >= 151,3) = 0;
subplot(2, 2, [1;3]);
scatter3(x_all,y_all,z_all,10,trialBlocks,'filled'); hold on;
% plot3(x_all(trial_ds<51),y_all(trial_ds<51),z_all(trial_ds<51),'b.');hold on; 
% plot3(x_all(trial_ds>50 & trial_ds<101),y_all(trial_ds>50 & trial_ds<101),z_all(trial_ds>50 & trial_ds<101),'k.');
% plot3(x_all(trial_ds>100 & trial_ds<151),y_all(trial_ds>100 & trial_ds<151),z_all(trial_ds>100& trial_ds<151),'r.');
% plot3(x_all(trial_ds>150),y_all(trial_ds>150),z_all(trial_ds>150),'g.'); 
xlabel('PC1'); 
ylabel('PC2');
zlabel('PC3');
t=title(filename);
set(t,'Interpreter','none')
% legend('Baseline','ControlInjx','Ketamine:1-50','Ketamine:51-200');
hold off;

subplot(2, 2, 4);
posx_ds = downsample(posx, ds_factor); 
scatter3(x_all,y_all,z_all,10,posx_ds,'filled'); hold on;
colorbar;

% plot(x_all(trial_ds>50 & trial_ds<101),y_all(trial_ds>50 & trial_ds<101),'k-');
% plot(x_all(trial_ds>100 & trial_ds<151),y_all(trial_ds>100 & trial_ds<151),'r-');
% plot(x_all(trial_ds>150),y_all(trial_ds>150),'g-'); 
lgd = legend('Baseline','ControlInjx','Ketamine:1-50','Ketamine:51-200');
lgd.FontSize = 18;
ch = findobj(get(lgd,'children'), 'type', 'line'); %// children of legend of type line
set(ch, 'Markersize', 12); %// set value as desired
xlabel('PC1'); 
ylabel('PC2');
hold off;

% subplot(2, 2, 4);
% plot(x_all(trial_ds<51),y_all(trial_ds<51),'b-'); hold on;
% plot(x_all(trial_ds>50 & trial_ds<101),y_all(trial_ds>50 & trial_ds<101),'k-');
% plot(x_all(trial_ds>100 & trial_ds<151),y_all(trial_ds>100 & trial_ds<151),'r-');
% plot(x_all(trial_ds>150),y_all(trial_ds>150),'g-'); 
% lgd = legend('Baseline','ControlInjx','Ketamine:1-50','Ketamine:51-200');
% lgd.FontSize = 18;
% ch = findobj(get(lgd,'children'), 'type', 'line'); %// children of legend of type line
% set(ch, 'Markersize', 12); %// set value as desired
% xlabel('PC1'); 
% ylabel('PC2');
% hold off;
%%
figure(5); hold on;
plot(coeff(:, 1))
plot(coeff(:, 2))
plot(coeff(:, 3))
%%
saveName = fullfile(imgDir, strcat(filename,'_PCA.jpg'));
saveas(gcf,saveName, 'jpg');
end