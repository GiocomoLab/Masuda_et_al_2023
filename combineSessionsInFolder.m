

indvSeshDirPath = '/Users/KeiMasuda/Desktop/indvSesh';
indvSeshes =  dir(fullfile(indvSeshDirPath,'*.mat'));
sessionMetaData = readtable('/Users/KeiMasuda/Desktop/fkm_analysis/SessionList.xlsx');

n = 50;
n = 69;
for n = 1:size(sessionMetaData,1)
    checkName = string(sessionMetaData{n,1});
    idx = contains({indvSeshes.name},checkName);
    seshesToCombine = indvSeshes(idx);
    try
        data_dir = seshesToCombine(1).folder;
        seshIDs =  extractfield(seshesToCombine,'name');
        
        % sort sheshIDs alphabetically (assumes alphabetical naming of
        % sessions is equivalent to real order)
        [~, neworder] = sort(lower(seshIDs));
        seshIDs = seshIDs(:,neworder);
        
        seshIDparts = cellfun(@returnStrippedSeshID,seshIDs,'UniformOutput',false);
        session_name = strcat(checkName, '_',strjoin(seshIDparts,'+'));
        unitySessions = cellfun(@removeMat, seshIDs,'UniformOutput',false);
        
        %Stitch file tother 
        stitchSynchedNPdata(data_dir,session_name, unitySessions); %unitySessions is a cell array of session names; can be a cell array of one name
        fprintf('Finished %d/%d\n',n,size(sessionMetaData,1));
    catch
       fprintf('Could not stitch %s\n',checkName) 
    end

end