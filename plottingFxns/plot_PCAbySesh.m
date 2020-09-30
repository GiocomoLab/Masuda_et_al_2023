function plot_PCAbySesh(cells,savefig)
%% Plot PCA by sessions
seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));

for i = 1:numel(seshes)
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(cells,seshIndx);
    try
        plot_PCA(seshCells,savefig)
        fprintf('Finished sesssion: %d/%d\n',i,numel(seshes))
        if ~savefig
            pause
        end
    catch
        fprintf('FAILED sesssion: %d/%d\n',i,numel(seshes))
    end
end

