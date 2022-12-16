addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat');
%%
allCells_fr = [];
allCells_trialDS = [];
for n = 1:numel(sessions)
    try
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        load(matPath);
        if size(all_fr,2) ~= 300
            fprintf('Not 300 trials for this session\n');
            continue
        end
        allCells_fr = vertcat(allCells_fr, all_fr);
        trial_ds = downsample(trials, ds_factor); 
        allCells_trialDS = vertcat(allCells_trialDS, trial_ds);
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end
%%
allCells_fr_pca = reshape(allCells_fr, size(allCells_fr,1),[]);
[coeff, ~, ~, ~, explained] = pca(allCells_fr_pca);
clf;
figure(1);
set(0,'defaultAxesFontSize',25)
set(gcf,'Position',[100 100 2000 1000]); 
subplot(2, 2, 2);
plot(cumsum(explained(1:10)));
title('% Variance Explained');
xlabel('Principle Components');
ylabel('% Variance Explained');

% Project onto the first three PCA components and plot as an animated
% line
firstThreeComponentsAll = coeff(:, [1:3]);

projectedDataAll = allCells_fr_pca * firstThreeComponentsAll;

x_all=projectedDataAll(:,1);
y_all=projectedDataAll(:,2);
z_all=projectedDataAll(:,3);


subplot(2, 2, [1;3]);
plot3(x_all,y_all,z_all,'r.');hold on; 
xlabel('PC1'); 
ylabel('PC2');
zlabel('PC3');
t=title(filename);
set(t,'Interpreter','none')
% legend('Baseline','ControlInjx','Ketamine:1-50','Ketamine:51-200');
hold off;

subplot(2, 2, 4);
plot(x_all,y_all,'b.'); hold on;
lgd = legend('Baseline','ControlInjx','Ketamine:1-50','Ketamine:51-200');
lgd.FontSize = 18;
ch = findobj(get(lgd,'children'), 'type', 'line'); %// children of legend of type line
set(ch, 'Markersize', 12); %// set value as desired
xlabel('PC1'); 
ylabel('PC2');
hold off;