% Before this, run masterScript to generate allCells
% Filter for WT mec cells during ketamine sessions
seshIndx = ismember(allCells.metadata(:,8),'ketamine');
ketamineCells = filterAllCellsStruct(allCells,seshIndx);
seshIndx = ismember(ketamineCells.metadata(:,4),'WT');
wt_ket_Cells = filterAllCellsStruct(ketamineCells,seshIndx);
fprintf('done filtering for WT-ket cells\n');

%% Session list & injection times synced to preprocessed spike times & position
wt_ket_sessions = cell(30,4);
session = '';
s = 1;
sess_idx = 0;
while s<=length(wt_ket_Cells.metadata(:,1))
    while strcmp(session, wt_ket_Cells.metadata(s,1))
        s = s+1;
    end
    session = wt_ket_Cells.metadata(s,1);
    sess_idx = sess_idx+1;
    substr = regexp(wt_ket_Cells.metadata(s,1),'_','split');
    wt_ket_sessions{sess_idx,1} = substr{1}{1};
    substr = regexp(wt_ket_Cells.metadata(s,1),'^.*_\d{6}','match');
    wt_ket_sessions{sess_idx,2} = char(substr{1});
    wt_ket_sessions{sess_idx,3} = wt_ket_Cells.posT(s).post(find(wt_ket_Cells.trial(s).trial==51, 1, 'first'));
    wt_ket_sessions{sess_idx,4} = wt_ket_Cells.posT(s).post(find(wt_ket_Cells.trial(s).trial==101, 1, 'first'));
end
save('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\wt_ket_sessions.mat', 'wt_ket_sessions')

%% add times synced to raw data
load('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\wt_ket_sessions.mat', 'wt_ket_sessions')
for s = 1:length(wt_ket_sessions)
	data_dir = ['\\oak-smb-giocomo.stanford.edu\groups\giocomo\export\data\Projects\JohnKei_NPH3\',...
        wt_ket_sessions{s,1},'\',wt_ket_sessions{s,2},'_keicontrasttrack_ketamine1_g0'];
	base_name = erase(wt_ket_sessions{s,2},[wt_ket_sessions{s,1},'_']);
    session_name = [base_name,'_keicontrasttrack_baseline1'];
    if ~isfile(fullfile(data_dir,strcat(session_name,'_position.txt')))
        base_name = [wt_ket_sessions{s,1},'_',base_name];
        session_name = [base_name,'_keicontrasttrack_baseline1'];
    end
    [wt_ket_sessions{s,5}, wt_ket_sessions{s,6}] = get_vr_start_end(data_dir, session_name);
    session_name = [base_name,'_keicontrasttrack_controlinjx1'];
    [wt_ket_sessions{s,7}, wt_ket_sessions{s,8}] = get_vr_start_end(data_dir, session_name);
    session_name = [base_name,'_keicontrasttrack_ketamine1'];
    [wt_ket_sessions{s,9}, wt_ket_sessions{s,10}] = get_vr_start_end(data_dir, session_name);
end
save('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\wt_ket_sessions.mat', 'wt_ket_sessions')