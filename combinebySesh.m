function combinebySesh(cells,foldername)
% Plot rasters from a single session as a combined image
% Inputs: allCellsStruct with session metadata in the first column
seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));
size_vert = 1042; % changed from 1042 to 346
size_horiz = 333; % changed from 333 to 667 for repeating tracks
numrow = 6; % number of rows in final image
for i = 1:numel(seshes)
    session_name = seshes{i};
    try
        
        combineTopDownRASTERS(session_name,foldername, size_vert, size_horiz, numrow)

        fprintf('Finished sesssion: %d/%d\n',i,numel(seshes))
    catch
        fprintf('FAILED sesssion: %d/%d\n',i,numel(seshes))
    end
end

