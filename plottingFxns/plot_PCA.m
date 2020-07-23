function plot_PCA(cells,savefig)
nCells = size(cells.spatialFR2,1);

all_fr = [];
ds_factor = 100;

trial = squeeze(cell2mat(struct2cell(cells.trial)));
trial = trial(:,1);
trial_ds = downsample(trial, ds_factor); 



% gaussian filter for smoothing
smoothSigma = 10;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

% Calculating smoothed firing rate
for i = 1:nCells
    cellFR = cells.FRtime(i).FRtime;
    % smooth firing rate
    fr_smoothed = conv(repmat(cellFR,1,3),gauss_filter,'same');
    fr = fr_smoothed(numel(cellFR)+1:numel(cellFR)*2);
    %
    all_fr(i, :) = fr;
end 

all_fr = all_fr';
all_fr_ds = downsample(all_fr, ds_factor);


%%
zFR = zscore(all_fr_ds);
[coeff, score, ~, ~, explained] = pca(zFR);
clf;


x_all=score(:,1);
y_all=score(:,2);
z_all=score(:,3);
%%
figure(1);
set(gcf,'Position',[100 100 1000 1000]); 


%
subplot(2,2,1)
% scatter3(x_all,y_all,z_all,10,trialBlocks,'filled'); hold on;
plot3(x_all(trial_ds<51),y_all(trial_ds<51),z_all(trial_ds<51),'k.');hold on; 
plot3(x_all(trial_ds>50 & trial_ds<101),y_all(trial_ds>50 & trial_ds<101),z_all(trial_ds>50 & trial_ds<101),'m.');
plot3(x_all(trial_ds>100 & trial_ds<151),y_all(trial_ds>100 & trial_ds<151),z_all(trial_ds>100& trial_ds<151),'b.');
plot3(x_all(trial_ds>150),y_all(trial_ds>150),z_all(trial_ds>150),'g.'); 
xlabel('PC1'); 
ylabel('PC2');
zlabel('PC3');
title(sprintf('%s, %s, %s, %s', cells.metadata{1,2},cells.metadata{1,3},cells.metadata{1,4},cells.metadata{1,8}));

set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');

subplot(2,2,2)
set(0,'defaultAxesFontSize',25)
figure(1)
plot(cumsum(explained(1:10)));
title('% Variance Explained');
xlabel('Principal Components');
ylabel('% Variance Explained');
box off;
axis square;

subplot(2,2,3)
plot(x_all(trial_ds<51),y_all(trial_ds<51),'k.'); hold on;
plot(x_all(trial_ds>50 & trial_ds<101),y_all(trial_ds>50 & trial_ds<101),'m.');
plot(x_all(trial_ds>100 & trial_ds<151),y_all(trial_ds>100 & trial_ds<151),'b.');
plot(x_all(trial_ds>150),y_all(trial_ds>150),'g.'); 
xlabel('PC1'); 
ylabel('PC2');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');

subplot(2,2,4)
plot(z_all(trial_ds<51),y_all(trial_ds<51),'k.'); hold on;
plot(z_all(trial_ds>50 & trial_ds<101),y_all(trial_ds>50 & trial_ds<101),'m.');
plot(z_all(trial_ds>100 & trial_ds<151),y_all(trial_ds>100 & trial_ds<151),'b.');
plot(z_all(trial_ds>150),y_all(trial_ds>150),'g.'); 
xlabel('PC3'); 
ylabel('PC2');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');

hold off;

if savefig
   imgDir = '/Users/KeiMasuda/Desktop/fkm_analysis/pca';
   saveName = fullfile(imgDir, strcat(sprintf('%s_%s_%s_%s', cells.metadata{1,2},cells.metadata{1,3},cells.metadata{1,4},cells.metadata{1,8}),'pca.jpg'));
   saveas(gcf,saveName, 'jpg'); 
end

end