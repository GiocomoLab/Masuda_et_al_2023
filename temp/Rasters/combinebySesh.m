function combinebySesh(cells,foldername, size_vert, size_horiz)
% Plot rasters from a single session as a combined image
% Inputs: allCellsStruct with session metadata in the first column
seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));
if isempty(size_vert)
    size_vert = 1042; % changed from 1042 to 346
end
if isempty(size_horiz)
    size_horiz = 333; % changed from 333 to 667 for repeating tracks
end

for i = 1:numel(seshes)
    session_name = seshes{i};
    try
        combineTopDownRASTERS(session_name,foldername, size_vert, size_horiz)

        fprintf('Finished sesssion: %d/%d\n',i,numel(seshes))
    catch
        fprintf('FAILED sesssion: %d/%d\n',i,numel(seshes))
    end
end

