
close all;
figure(12); hold on;
plot(nanmean(all_cellCorrScore_WT(:,1:300),1),'LineWidth',5,'DisplayName','WT')
plot(nanmean(all_cellCorrScore_KO(:,1:300),1),'LineWidth',5,'DisplayName','HCN1 KO')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Average Correlation Score Curve WT vs KO'))
xlabel('Trial')
ylabel('Corrleation compared to Baseline Template (rho)')
legend;

addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
%%

avgControlInjxWT = nanmean(nanmean(all_cellCorrScore_WT(:,55:100),1),2);
avgRestingWT = nanmean(all_cellCorrScore_WT(:,101:290),1);
stabilityRatioWT = avgRestingWT/avgControlInjxWT;


avgControlInjxKO = nanmean(nanmean(all_cellCorrScore_KO(:,55:100),1),2);
avgRestingKO = nanmean(all_cellCorrScore_KO(:,101:290),1);
stabilityRatioKO = avgRestingKO/avgControlInjxKO;


close all;
figure(12); hold on;
plot(stabilityRatioWT,'LineWidth',5,'DisplayName','WT');
plot(stabilityRatioKO,'LineWidth',5,'DisplayName','HCN1 KO')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Correlation Score Curve Normalized by ControlInjx'))
xlabel('Trial')
ylabel('Corrleation compared to Baseline Template (rho)')
legend;