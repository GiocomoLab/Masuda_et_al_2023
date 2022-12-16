function plot_correlationMatrixBySessionThenAvg(image_save_dir, cells,save_figs)
close all;
% image_save_dir = '/Users/fkmasuda/Desktop/fkm_analysis/correlation_matrices'
seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));

corrmatrices = nan(numel(seshes),290,290);
for i = 1:numel(seshes)
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    sesh_cells = filterAllCellsStruct(cells,seshIndx);
    corrMatrix = calc_correlationMatrix(sesh_cells);
    corrmatrices(i,:,:) = corrMatrix;
    
    clc;
    h = figure(1);
    imagesc(corrMatrix);
%     imagesc(corrMatrix,[0.5 0.9]);
    goodFigPrefs;
    colormap('bone');
    colorbar;
    title_name = sprintf('Sesh %i',i);
    title(title_name);
    if save_figs
            saveas(h,fullfile(image_save_dir,sprintf('%s.png',title_name)),'png');
    else
        pause
    end
end
avg_corrmatrices = squeeze(nanmean(corrmatrices,1));

figure()
imagesc(avg_corrmatrices);
goodFigPrefs;
colormap('bone');
colorbar;
title('Average Session Spatial Correlation')