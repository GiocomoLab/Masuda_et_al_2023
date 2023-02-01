filenames = dir('Z:\giocomo\export\data\Projects\JohnKei_NPH3\**\*.mat');
%%
% sessionList = {'G1_190817','G1_190818','G2_190712','G2_190714','G2_190715','G2_190716','G3_190703','G3_190704','G3_190705','G3_190708','G4_190619','G4_190620','G4_190621','G4_190623','G4_190624','G5_190703','G5_190704','G5_190705','G5_190708','HCN1_190618','HCN1_190619','HCN1_190621','HCN1_190623','HCNb2_191010','HCNb2_191012','HCNb2_191014','HCNb2_191015','HCNb2_191016','HCNb4_191010','HCNb4_191012','HCNd1_190808','HCNd1_190809','HCNd1_190812','HCNd1_190813','HCNd2_190820','HCNd2_190821','HCNe1_190824','HCNe1_190827','HCNe1_190828','HCNe1_190830','HCNe1_190831','HCNe2_190912','HCNe2_190913','HCNe2_190914','HCNe3_190929','HCNe3_190930','HCNe3_191001','HCNe3_191002','HCNe3_191003','npI1_190418','HCNf1_191214','HCNf1_191216','HCNf1_191219','HCNf1_191220','HCNf2_191214','HCNf2_191215','HCNf2_191218','HCNf2_191220','HCNg2_191210','HCNg2_191212','HCNg2_191215','HCNh3_200207','HCNh3_200212','HCNh3_200213','HCNh3_200214','HCNi1_200212','HCNi1_200218','HCNi1_200219','HCNi1_200220','HCNi2_200130','HCNi2_200131','HCNi2_200202','HCNi2_200203','HCNi2_200204','HCNi2_200205','HCNj2_200130','HCNj2_200131','HCNj2_200202','HCNj2_200203','HCNj2_200204'};
sessionList = {'G3_190703','G3_190704','G4_190620','HCNe2_190914'};
sessionList = {'npI1_190418'}
% sessionList = {'npI1_190418','HCNi1_200220'};
idx = contains({filenames.folder},sessionList);
fltrFilenames = filenames(idx);
 
removeList = {'+','rez','neuropix'};
rmidx = contains({fltrFilenames.name},removeList);
fltrFilenames2 = fltrFilenames(~rmidx);
%%
 
for i = 1:numel(fltrFilenames2)
    tracklength = 400;
    save_fig_bool = true;
    save_data_dir = 'Z:\giocomo\export\data\Projects\JohnKei_NPH3\data_giocomo\individualUnitySessionData';
 
    % location of data
    data_dir = fltrFilenames2(i).folder;
    [~,name,~] = fileparts(fltrFilenames2(i).name);
    session_name = name(4:end);
    if isfile(fullfile(data_dir,strcat(session_name,'_position.txt')))
        session_name = name(4:end);
    elseif isfile(fullfile(data_dir,strcat(name(5:end),'_position.txt')))
        session_name = name(5:end);
    elseif isfile(fullfile(data_dir,strcat(name(6:end),'_position.txt')))
        session_name = name(6:end);
    elseif isfile(fullfile(data_dir,strcat(name(7:end),'_position.txt')))
        session_name = name(7:end);
    elseif isfile(fullfile(data_dir,strcat(name(8:end),'_position.txt')))
        session_name = name(8:end);
    elseif isfile(fullfile(data_dir,strcat(name(9:end),'_position.txt')))
        session_name = name(9:end);
    elseif isfile(fullfile(data_dir,strcat(name(10:end),'_position.txt')))
        session_name = name(10:end);
    elseif isfile(fullfile(data_dir,strcat(name,'_position.txt')))
        session_name = name;
    else
        fprintf('Could not find position file for %s\n',session_name);
        continue
    end
    try
        sync_vr_to_np_driftCorrect(data_dir, save_data_dir, session_name, tracklength, save_fig_bool)
    catch
        fprintf('Could not sync %s\n',session_name)
        fclose('all')
    end
end
