function plot_avgFRbyMouse(allCells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Firing Rate over Trials by Mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seshes = unique(cellfun(@num2str,allCells.metadata(:,2),'uni',0));

figure(); hold on;
for i = 1:numel(seshes)
    
    seshIndx = ismember(allCells.metadata(:,2),seshes{i});
    cellsFR = allCells.spatialFR10(seshIndx,:,:);
    plot(smooth(nanmean(nanmean(cellsFR,3),1),10),'-k');
%     plot_lineWithSEM(nanmean(cellsFR,3),[])
end
plot(smooth(nanmean(nanmean(allCells.spatialFR10,3),1)),'-r','LineWidth',5)
title('Average FR by Mouse');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
ylabel('Hz')
xlabel('Trial Number')
% plot_lineWithSEM(nanmean(allCells.spatialFR10,3),[])
end

