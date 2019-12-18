function plot_HistfitKetCorrEffectScore(allCells, cellIndx, filter)
% Plot Distribution of Ketamine Correlation Scores
% Input: allCells struct, cell index to filter which cells to plot, name of
% filter

figure();
des = allCells.drugEffectScores(cellIndx,:);

hold on;
histfit(allCells.drugEffectScores((des(:,4) < 2),4),50)
%         histfit(allCells.drugEffectScores(:,2),50, 'kernel')
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
title(sprintf('Distribution of Ketamine Correlation Effect Scores(%s)',filter))
xlabel('Pre vs Post Ketamine Correlation Effect Score')
ylabel('Number of Cells')
vline(0,'r')

end