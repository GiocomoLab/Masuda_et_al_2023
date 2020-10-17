function plot_correlationMatrixBySessionThenAvg(cells)

seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));

corrmatrices = nan(numel(seshes),290,290);
for i = 1:numel(seshes)
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    sesh_cells = filterAllCellsStruct(cells,seshIndx);
    corrmatrices(i,:,:) = calc_correlationMatrix(sesh_cells);
    
%     
%     figure()
%     imagesc(avg_corrmatrices,[0.5 0.9]);
%     goodFigPrefs;
%     colormap('bone');
%     colorbar;
%     title(sprintf('Sesh %i',i))
end
avg_corrmatrices = squeeze(nanmean(corrmatrices,1));

figure()
imagesc(avg_corrmatrices,[0.5 0.9]);
goodFigPrefs;
colormap('bone');
colorbar;
title('Average Session Spatial Correlation')