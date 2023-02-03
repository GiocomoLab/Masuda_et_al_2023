function plot_BehaviorbySesh(cells,savefigTF)
%% Plot Lick Behavior by sessions
% plots lick accuracy over trials

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
%% Plot Stats of Mean lick accuracies of the first 50 trials of each session
figure(2)
x = horzcat(...
    repmat({'Baseline'}, [1 50]),...
    repmat({'Control'}, [1 50]),...
    repmat({'Ketamine'}, [1 50])...
    );

y = nanmean(allLickAccuracy(:,1:150),1);

customColorMap = [ 0.5 0.5 0.5 %grey
    0.8 0.2 0.8 %magenta
    0 0.8 0.2]; %green

clear g;
g=gramm('x',x,'y',y,'color',color);
g.stat_violin('fill', 'transparent');
% g.stat_boxplot('notch','true')
% g.stat_summary('geom',{'edge_bar','black_errorbar'},'type','ci','setylim','true');
g.geom_jitter('width',0.3,'height',0,'dodge',1,'alpha',1);
g.set_color_options('chroma',0,'map',customColorMap);
g.set_names('x','','y', 'Lick Accuracy');
g.set_order_options('x',0);
g.draw();

figure(3)
[p,t,stats] = anova1(y,x,'off');
[results,~] = multcompare(stats,'CType','bonferroni');




