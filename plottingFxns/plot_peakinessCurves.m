function plot_peakinessCurves(allCells)
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
metadata = allCells.metadata; %session_name, animalName, sessionDate, genotype, gender, ketamine_day
% Genotype Stairs
close all;
figure();
%%
genotype = 'WT';
category = 'peakiness';

% subplot(3,3,[1,4,7])
% WTindx = ismember(metadata(:,4),genotype);
% plot_sortedMatrix(allCells.(category),WTindx,'descend')
% colorbar;
% title('Peakiness');

subplot(3,3,[1, 2,3]);
plot_lineWithSEM(allCells.(category),[])
title('Crest Factor');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');

category = 'bitsPerSpikeCurve';
subplot(3,3,[4,5,6]);
plot_lineWithSEM(allCells.(category),[])
title('Bits Per Spike');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');

category = 'stabilityScoreCurve';
subplot(3,3,[7,8,9]);
plot_lineWithSEM(allCells.(category),[])
title('Stability');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
end