function plot_HistfitKetCorrEffectScore(allCells, filter)
% Plot Distribution of Ketamine Correlation Scores
% Input: allCells struct, cell index to filter which cells to plot, name of
% filter

figure();
des = allCells.drugEffectScores;

hold on;
h = histfit(allCells.drugEffectScores((des(:,4) < 2),4),50);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Distribution of Ketamine Correlation Effect Scores(%s)',filter))
xlabel('Normalized Correlation to Baseline Template (rho)')
ylabel('Number of Cells')
vline(0,'r')
h(1).FaceColor = [0.8 0.8 0.8];
h(2).Color = [.1 .1 .1];
end