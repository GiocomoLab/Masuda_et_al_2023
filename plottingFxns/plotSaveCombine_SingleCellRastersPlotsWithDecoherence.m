function plotSaveCombine_SingleCellRastersPlotsWithDecoherence(cells,dch,image_save_dir,save_figs)
% Plot Single Cell Raster plots with Decoherence period 
% Highlighted with trial-based Dch Index
% Input: allCells struct, dch struct

% Find unique sessions in cells struct
seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));
for i = 1:numel(seshes)
    close all;
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(cells,seshIndx);
    
    plotDchRaster(dch.decoherenceIdx,seshCells,i, save_figs,image_save_dir)
end
% Combine single cell rasters into large session pngs 
combinebySesh(fltrCells,image_save_dir)

end