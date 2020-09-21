function plot_avgFRbyMouse(allCells,titleStr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Firing Rate over Trials by Mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seshes = unique(cellfun(@num2str,allCells.metadata(:,2),'uni',0));

figure(); hold on;
for i = 1:numel(seshes)
    
    seshIndx = ismember(allCells.metadata(:,2),seshes{i});
    cellsFR = allCells.spatialFRsmooth(seshIndx,:,:);
    plot(smooth(nanmean(nanmean(cellsFR,3),1),10),'-k');
%     plot_lineWithSEM(nanmean(cellsFR,3),[])
end
plot(smooth(nanmean(nanmean(allCells.spatialFRsmooth,3),1)),'-r','LineWidth',5)
title(titleStr);
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
% axis square;
set(gca,'FontSize',20);
set(gca,'FontName','Helvetica');
ylabel('Hz')
xlabel('Trial Number')
% plot_lineWithSEM(nanmean(allCells.spatialFR10,3),[])

vline(50,'m','Control Injx')
vline(100,'g','Ketamine Injx')
end

