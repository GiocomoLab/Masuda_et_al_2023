function plot_avgFRbyMouse(cells,titleStr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Firing Rate over Trials by Mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seshes = unique(cellfun(@num2str,cells.metadata(:,2),'uni',0));

figure(); hold on;

for i = 1:numel(seshes)
    
    seshIndx = ismember(cells.metadata(:,2),seshes{i});
    cellsFR = cells.spatialFRsmooth(seshIndx,:,:);
    
    avgFRline = smooth(nanmean(nanmean(cellsFR,3),1),10);
    baselineFR = nanmean(nanmean(nanmean(cellsFR(:,1:50,:),3),1),2);
    normalizedFRline = avgFRline-baselineFR;
    plot(normalizedFRline(1:290),'-k');
    normalizedLines(i,:,:) = normalizedFRline(1:290);
%     plot_lineWithSEM(nanmean(cellsFR,3),[])
end
% plot(smooth(nanmean(nanmean(allCells.spatialFRsmooth,3),1)),'-r','LineWidth',5)

plot_lineWithSEM(nanmean(normalizedLines,3),[])
plot(nanmean(normalizedLines,1),'-r','LineWidth',4)

% title(titleStr);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
% box off;
% axis square;
set(gca,'FontSize',40);
set(gca,'FontName','Helvetica');
ylabel('Mean FR Change(Hz)')
xlabel('Trial Number')
% plot_lineWithSEM(nanmean(allCells.spatialFR10,3),[])

vline(50,'m','Control Injx')
vline(100,'g','Ketamine Injx')
end

