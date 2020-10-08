function plot_peakinessCurves(cells)
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
metadata = cells.metadata; %session_name, animalName, sessionDate, genotype, gender, ketamine_day
% Genotype Stairs
close all;
figure();
category = 'peakiness';
trialMax = 290;

% subplot(3,3,[1,4,7])
% WTindx = ismember(metadata(:,4),genotype);
% plot_sortedMatrix(cells.(category),WTindx,'descend')
% colorbar;
% title('Peakiness');

subplot(3,3,[1,2,3]);
data = cells.(category);
plot_lineWithSEM(data(:,1:trialMax),[])
title('Crest Factor');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');

category = 'bitsPerSpikeCurve';
subplot(3,3,[4,5,6]);
data = cells.(category);
plot_lineWithSEM(data(:,1:trialMax),[])
title('Bits Per Spike');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');

category = 'stabilityScoreCurve';
subplot(3,3,[7,8,9]);
data = cells.(category);
plot_lineWithSEM(data(:,1:trialMax),[])
title('Stability');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
end