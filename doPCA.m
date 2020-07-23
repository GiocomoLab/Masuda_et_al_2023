function doPCA(matPath)

addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));
%% 
% matPath = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190620_keicontrasttrack_ketamine1_g0/G4_190620_keicontrasttrack_baseline+cntrlinjx+ketamine';
% load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/HCN1_190619_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');
% load('/Users/KeiMasuda/Dropbox/3_GiocomoLab/CodeGiocomoLab/githubRepos/JohnKeiNPAnalysis/logisticRegression/G3_190708_keicontrasttrack_baseline+cntrlinjx+ketamine.mat');

load(matPath);
%%
[cells_to_plot, ~,~] = spPreProcess(sp);
nCells =  numel(cells_to_plot);
%% Calculating firing rate

% calculate firing rate by time
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

%% PCA 
imgDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/img';
[~,filename] = fileparts(matPath);
%plot(mean(all_fr)
zFR = zscore(all_fr);
[coeff, ~, ~, ~, explained] = pca(zFR);
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

projectedDataAll = zFR * firstThreeComponentsAll;

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
% 
% subplot(2, 2, 4);
% posx_ds = downsample(posx, ds_factor); 
% scatter3(x_all,y_all,z_all,10,posx_ds,'filled'); hold on;
% colorbar;
% lgd = legend('Baseline','ControlInjx','Ketamine:1-50','Ketamine:51-200');
% lgd.FontSize = 18;
% ch = findobj(get(lgd,'children'), 'type', 'line'); %// children of legend of type line
% set(ch, 'Markersize', 12); %// set value as desired
% xlabel('PC1'); 
% ylabel('PC2');
% hold off;

subplot(2, 2, 4);
plot(x_all(trial_ds<51),y_all(trial_ds<51),'b.'); hold on;
plot(x_all(trial_ds>50 & trial_ds<101),y_all(trial_ds>50 & trial_ds<101),'k.');
plot(x_all(trial_ds>100 & trial_ds<151),y_all(trial_ds>100 & trial_ds<151),'r.');
plot(x_all(trial_ds>150),y_all(trial_ds>150),'g.'); 
lgd = legend('Baseline','ControlInjx','Ketamine:1-50','Ketamine:51-200');
lgd.FontSize = 18;

ch = findobj(get(lgd,'children'), 'type', 'line'); %// children of legend of type line
set(ch, 'Markersize', 12); %// set value as desired
xlabel('PC1'); 
ylabel('PC2');
axis('square');
box('off');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);mes
set(gca,'layer','bottom');

hold off;
%%
saveName = fullfile(imgDir, strcat(filename,'_PCA.jpg'));
saveas(gcf,saveName, 'jpg');
end