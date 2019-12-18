function plot_stabilityScore(allCells, cellsIndx,title)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stability Score Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(); clf;
clear g;
x = 1:290;
if isempty(cellsIndx)
    y = allCells.stabilityScoreCurve(:,x);
else
    y = allCells.stabilityScoreCurve(cellsIndx,x);
end
g(1,1) = gramm('x',x,'y',y); 
g(1,1).stat_summary('type','sem','setylim',true);

g(1,1).set_title(sprintf('Ketamine-induced Stability over Time: %s', title), 'FontSize', 40);
g(1,1).set_names('x','Trial','y','Stability');
g(1,1).set_text_options('base_size',20);
g.draw()
axis square;
set(gcf,'Position',[100 100 1000 1000])
set(gca,'TickDir','out');
end