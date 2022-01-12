function findAttractorDynamicsbySession(cells)

seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));
slopeMatrix = nan(numel(seshes),8);
for i = 1:numel(seshes)
    
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(cells,seshIndx);
    sf.name = seshCells.metadata{1,2};
    sf.genotype = seshCells.metadata{1,4};
    sf.sessionDate = seshCells.metadata{1,3};
    sf.image_save_dir = '../fkm_analysis/attractorNetworkCorrelations';
    sf.seshNum=i;
    try
        slopeMatrix(i,:) = findAttractorDynamics(seshCells,true,sf);
    catch
        fprintf(strcat("could not do"+num2str(i)+"\n"))
    end
end
%%
close all; clear g;
figure(); hold on;
linearizedSlopeMatrix = squeeze(reshape(slopeMatrix,[],1));

y = linearizedSlopeMatrix;

x = vertcat(...
    repmat({'baselineA vs baselineB'}, size(slopeMatrix,1),1),...
    repmat({'baseline vs cntrl'}, size(slopeMatrix,1),1),...
    repmat({'baseline vs acuteKet'}, size(slopeMatrix,1),1),...
    repmat({'cntrl vs acuteKet'}, size(slopeMatrix,1),1),...
    repmat({'baseline vs acuteKet'}, size(slopeMatrix,1),1),...
    repmat({'acuteKet vs lateKet'}, size(slopeMatrix,1),1),...
    repmat({'lateKet vs gainChange'}, size(slopeMatrix,1),1),...
    repmat({'cntrl vs lateKet'}, size(slopeMatrix,1),1)...
    );

clf;
figure(1);
clear g;
g=gramm('x',x,'y',y);
g.stat_violin('fill', 'transparent');
% g.stat_boxplot('notch','true')
% g.stat_summary('geom',{'edge_bar','black_errorbar'},'type','ci','setylim','true');
g.geom_jitter('width',0.6,'height',0,'dodge',1,'alpha',1);
g.set_color_options('chroma',0);
g.set_title('Correlations between Stable Cell Paris');
g.set_names('x','','y', 'Correlation of FR Stability between Different Conditions');
% x={'baselineA vs baselineB','baseline vs cntrl','baseline vs acuteKet','cntrl vs acuteKet','cntrl vs lateKet','baseline vs acuteKet','acuteKet vs lateKet','lateKet vs gainChange'};
g.set_order_options('x',0);
g.draw();
