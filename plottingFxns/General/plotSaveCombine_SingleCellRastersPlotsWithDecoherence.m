function plotSaveCombine_SingleCellRastersPlotsWithDecoherence(cells,image_save_dir,save_figs)
% Plot Single Cell Raster plots with Decoherence period 
% Highlighted with trial-based Dch Index
% Input: allCells struct, dch struct

% Find unique sessions in cells struct
seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));
for i = 1:numel(seshes)
    close all;
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(cells,seshIndx);
    
    plotDchRaster(seshCells.dch(1).dch.decoherenceIdx,seshCells, save_figs,image_save_dir)
end
% Combine single cell rasters into large session pngs
size_vert = 1042; % changed from 1042 to 346
size_horiz = 333; % changed from 333 to 667 for repeating tracks
combinebySesh(cells,image_save_dir,size_vert,size_horiz)

end