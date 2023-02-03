function plot_spatialFR_PCA(cells,savefig)
nCells = size(cells.spatialFR2,1);

all_fr = cells.spatialFR2;
all_fr_trialAvg = squeeze(nanmean(all_fr)); 

%%
[idx,C] = kmeans(all_fr_trialAvg,2);
scatter3(all_fr_trialAvg(:,1), all_fr_trialAvg(:,2),all_fr_trialAvg(:,3), idx);   % plot three clusters with different colors
hold on;
plot3(C(:, 1), C(:, 2), C(:, 3), 'kx');   % plot centroids
%%
X = all_fr_trialAvg
eva = evalclusters(all_fr_trialAvg,'kmeans','CalinskiHarabasz','KList',1:10)
plot(eva)

Xred = nnmf(all_fr_trialAvg,2);

gscatter(Xred(:,1),Xred(:,2),{'1','2'})
xlabel('Column 1')
ylabel('Column 2')
grid on
%%
imagesc(all_fr_trialAvg)
%%
pcaCoordTrialRange = 100:150;
zFR = zscore(all_fr_trialAvg(pcaCoordTrialRange,:));
[coeff, score, ~, ~, explained] = pca(zFR,'Centered',false);
clf;


% x_all=score(:,1);
% y_all=score(:,2);
% z_all=score(:,3);



firstThreeComponentsAll = coeff(:, [1:3]);

projectedDataAll = all_fr_trialAvg(pcaCoordTrialRange,:) * firstThreeComponentsAll;

x_all=projectedDataAll(:,1);
y_all=projectedDataAll(:,2);
z_all=projectedDataAll(:,3);
%
figure(1);
set(gcf,'Position',[100 100 1000 1000]); 
plot3(x_all,y_all,z_all)

%%
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