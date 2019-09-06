function sessions = filterSessions(sessions, filter)
% filter sessions by subgroups
% Input:
%   filter - filter keyword: e.g. mec, WT,KO
% Output:
%   sessions - filtered session list

if strcmp(filter, 'mec')
    remove = {'AA','B1','B3','E1','E2','F3','propofol','MK801','190809'}; %check for this in session name and remove
elseif strcmp(filter, 'WT')
    % Remove KO + strange sessions
    remove = {'AA','B1','B3','E1','E2','F3','propofol','MK801','190809','HCNd2','HCNe1','HCN1'};
elseif strcmp(filter, 'KO')
    % Remove WT + strange sessions
    remove = {'AA','B1','B3','E1','E2','F3','propofol','MK801','190809','G1','G2','G3','G4','G5','HCNd1','npI1'}; 
else
    fprintf('Bad filter key. Choose: mec, WT, KO');
end

for z= 1:numel(remove)
    idx = ~cellfun('isempty',strfind({sessions.name},remove{z}));
    sessions(idx) = [];
end
