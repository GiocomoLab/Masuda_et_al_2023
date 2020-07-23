function plot_BehaviorbySesh(cells,savefigTF)
%% Plot Lick Behavior by sessions
seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));
allLickAccuracy = nan(numel(seshes),300);
for i = 1:numel(seshes)
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(cells,seshIndx);
    try
        lickAccuracyByTrial = plot_SingleSessionBehavior(seshCells,savefigTF);
        if size(lickAccuracyByTrial,2) > 300
            allLickAccuracy(i,1:300) = lickAccuracyByTrial(1:300);
        else
            allLickAccuracy(i,:) = lickAccuracyByTrial(:);
        end
        fprintf('Finished sesssion: %d/%d\n',i,numel(seshes))
    catch
        fprintf('FAILED sesssion: %d/%d\n',i,numel(seshes))
    end
end
close all;
plot_lineWithSEM(allLickAccuracy,[])
% title('Lick Accuracy');
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
axis square;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
xlabel('Trials'); 
ylabel('Lick Accuracy');
