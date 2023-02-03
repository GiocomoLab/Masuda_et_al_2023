function plot_AttractorDynamicsbySession_figure(cells,save_figs,folderName)

folderName = "time_wt_cell_stable_Figure";

close all;
seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));
% slopeMatrix = nan(numel(seshes),4);
pearson_rhoMatrix = nan(numel(seshes),4);

for i = 1:numel(seshes)
    fprintf("session %i/%i\n",i,numel(seshes))
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(cells,seshIndx);
    sf.name = seshCells.metadata{1,2};
    sf.genotype = seshCells.metadata{1,4};
    sf.sessionDate = seshCells.metadata{1,3};
    sf.image_save_dir = strcat('../fkm_analysis/attractorNetworkCorrelations/',folderName);
    sf.seshNum=i;
    save_figs = false;
    try
        [~,pearson_rho] = findAttractorDynamics_figure(seshCells,save_figs,sf);
%         slopeMatrix(i,:) = slope;
        pearson_rhoMatrix(i,:) = pearson_rho;
    catch
        fprintf(strcat("could not do"+num2str(i)+"\n"))
    end
end
%%
close all; clear g;

%
k = figure(1);
set(k,'Position',[100 100 1200 400]);
time_only_rhoMatrix = pearson_rhoMatrix(:,1:4);
linearizedPearson_rhoMatrix = squeeze(reshape(time_only_rhoMatrix,[],1));
y = linearizedPearson_rhoMatrix;


rhoMatrixSize = size(pearson_rhoMatrix,1);


x = vertcat(...
	repmat({'baseline(time) vs cntrl(time)'}, rhoMatrixSize,1),...
	repmat({'baseline(time) vs acuteKet(time)'}, rhoMatrixSize,1),...
	repmat({'baseline(time) vs lateKet(time)'}, rhoMatrixSize,1),...
	repmat({'lateKet(time) vs gainChange(time)'}, rhoMatrixSize,1)...
    );

clear g;
g=gramm('x',x,'y',y);
g.stat_violin('fill', 'transparent');
% g.stat_boxplot('notch','true')
% g.stat_summary('geom',{'edge_bar','black_errorbar'},'type','ci','setylim','true');
g.geom_jitter('width',0.3,'height',0,'dodge',1,'alpha',1);
g.set_color_options('chroma',0);
% g.set_title('Pearson Correlations between Stable Cell Pairs');
g.set_names('x','','y', 'Rho');
% x={'baselineA vs baselineB','baseline vs cntrl','baseline vs acuteKet','cntrl vs acuteKet','cntrl vs lateKet','baseline vs acuteKet','acuteKet vs lateKet','lateKet vs gainChange'};
g.set_order_options('x',0);
g.draw();
if save_figs
    saveas(k,fullfile(sf.image_save_dir,folderName),'svg');
end
%%
figure(2)
[p,t,stats] = anova1(y,x,'off');
[results,~] = multcompare(stats,'CType','bonferroni');

% print(results)
