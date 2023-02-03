function plot_HistfitFRscore(allCells, filter)
% PLOT Histfit on FR score
% Input: allCells struct, cell index to filter which cells to plot, name of
% filter

figure();
des = allCells.drugEffectScores;

threshold = 15;
ketFRscore = des(des(:,3) < threshold & des(:,3) > -threshold,3);
h = histfit(ketFRscore,50);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Distribution of Ketamine FR Effect Scores(%s)',filter))
xlabel('FR Effect Score:Ratio of FR change between Drug & Control')
ylabel('Number of Cells')
vline(0,'r')
h(1).FaceColor = [0.6 0.6 0.6];
h(2).Color = [.1 .1 .1];

end