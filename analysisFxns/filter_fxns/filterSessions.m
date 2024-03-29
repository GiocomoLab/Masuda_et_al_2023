function [sessions,sessionMetaData] = filterSessions(sessions, sessionMetaData, filter)
% filter sessions and sessionMetaData by subgroups
% Input:
%   filter - filter keyword: e.g. mec, WT,KO
% Output:
%   sessions - filtered session list

if strcmp(filter, 'mec')
    remove = {'AA','B1','B3','E1','E2','F3','propofol','all','john','D2','C2','HCN1_190620','G2_190702','john','Baseline','D2','npI1_190418','HCNi1_200220','HCNd1_190813'}; %check for this in session name and remove
elseif strcmp(filter, 'WT')
    % Remove KO + strange sessions
    remove = {'AA','B1','B3','E1','E2','F3','propofol','MK801','HCNd2','HCNe1','HCNe3','HCN1','all','john','Baseline','D2','HCNb4'};
elseif strcmp(filter, 'KO')
    % Remove WT + strange sessions
    remove = {'AA','B1','B3','E1','E2','F3','propofol','MK801','G1','G2','G3','G4','G5','HCNd1','HCNe2','npI1','all', 'john','Baseline','D2','HCNb2'}; 
elseif strcmp(filter, 'MK801')
    % Remove WT + strange sessions
    keep = {'MK801'}; 
    for z= 1:numel(keep)
        idx = cellfun('isempty',strfind({sessions.name},keep{z}));
        sessions(idx) = [];
    end
    return
elseif strcmp(filter, 'control')
    % Remove WT + strange sessions
    keep = {'control'}; 
    for z= 1:numel(keep)
        idx = cellfun('isempty',strfind({sessions.name},keep{z}));
        sessions(idx) = [];
    end
    return
else
    fprintf('Bad filter key. Choose: mec, WT, KO');
end

for z= 1:numel(remove)
    idx = ~cellfun('isempty',strfind({sessions.name},remove{z}));
    sessions(idx) = [];
    % if sessionMetaData is present; filter that too
    if ~isempty(sessionMetaData)
        cellSMD = table2cell(sessionMetaData);
        idx = ~cellfun('isempty',strfind(sessionMetaData.Session,remove{z}));
        sessionMetaData(idx,:) = [];
    else
        sessionMetaData = [];
    end
end
