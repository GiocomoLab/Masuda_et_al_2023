%% Load data, get spikes from good cells

load('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/F3/F3_190620_johncontrasttrack_train1_g0/F3_190620_johncontrasttrack_train1_incl_earlyend.mat');
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
end_recording = max(post);
pre_earlyend = [-60:0.02:end_recording+60]'; 
all_fr = [];
ds_factor = 20;
smoothSigma = 10;

for i = 1:nCells
   
    % get spike times for cell i
    spike_t = sp.st(sp.clu==cells_to_plot(i));
    
    % calculate firing rate
    fr = calcSmoothedFR_Time(pre_earlyend, spike_t, ds_factor, smoothSigma); 
    all_fr(i, :) = fr;
        
end

all_fr = all_fr';

time_pre_VR = 60; % assumes pre and post-VR times are equivalent

pre_VR_frames = time_pre_VR/(0.02*ds_factor);

VR_fr = all_fr(pre_VR_frames:size(all_fr, 1)-pre_VR_frames, :); % firing rates from times mouse was shown VR
pre_fr = all_fr(1:pre_VR_frames, :);

%% PCA 

nComponents = 3;

[eigvecs_all, transformed_all, eigvals_all, ~, ~, mu_all] = pca(all_fr);
[eigvecs_vr, transformed_vr, eigvals_vr, ~, ~, mu_vr] = pca(VR_fr);
[eigvecs_pre, transformed_pre, eigvals_pre, ~, ~, mu_pre] = pca(pre_fr);

% Scree plot to visually determine how many principal components are
% needed to capture variance in the data

figure();
plot(eigvals_all);
title('Eigenvalues on all data')

figure();
plot(eigvals_vr);
title('Eigenvalues on VR data')

figure()
plot(eigvals_pre);
title('Eigenvalues on pre-VR data')

% use the elbow in the plot to determinee how many principal components are
% needed. 


% Project onto the first three PCA components and plot as an animated
% line

firstThreeComponentsAll = eigvecs_all(:, [1:3]);
firstThreeComponentsVR = eigvecs_vr(:, [1:3]);
firstThreeComponentsPre = eigvecs_pre(:, [1:3]);

projectedDataAll = all_fr * firstThreeComponentsAll;
projectedDataVR = all_fr * firstThreeComponentsVR;
projectedDataPre = all_fr * firstThreeComponentsPre;

x_all=projectedDataAll(:,1);
y_all=projectedDataAll(:,2);
z_all=projectedDataAll(:,3);

x_vr=projectedDataVR(:,1);
y_vr=projectedDataVR(:,2);
z_vr=projectedDataVR(:,3);

x_pre=projectedDataPre(:,1);
y_pre=projectedDataPre(:,2);
z_pre=projectedDataPre(:,3);

% Code for plotting animated line
animatedPCAearly = animatedline('Color', 'red');
animatedPCAtrials = animatedline('Color', 'blue');
animatedPCAend = animatedline('Color', 'magenta');

axis([min(x_all) max(x_all) min(y_all) max(y_all) min(z_all) max(z_all)]);
view(3);
hold on;
title('Activity Along Top 3 Principal Components Before, During, and After VR')

for i = 1:pre_VR_frames
    addpoints(animatedPCAearly, x_all(i), y_all(i), z_all(i));
    drawnow
end

for i = pre_VR_frames+pre_VR_frames+1:length(x)-pre_VR_frames-pre_VR_frames
    addpoints(animatedPCAtrials, x_all(i), y_all(i), z_all(i));
    drawnow
end

for i = length(x)-pre_VR_frames-pre_VR_frames:length(x)
    addpoints(animatedPCAend, x_all(i), y_all(i), z_all(i));
    drawnow
end

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




