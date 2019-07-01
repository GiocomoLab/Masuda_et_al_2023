function stitchSynchedNPdata(data_dir, session_name, unitySessions)

% Input:
data_dir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190625_keicontrasttrack_propofol1_g0';
session_name = 'G4_190625_keicontrasttrack_baseline+cntrlinjx+propofol';
unitySessions ={'G4_190625_keicontrasttrack_baseline1', 'G4_190625_keicontrasttrack_controlinjx1','G4_190625_keicontrasttrack_propofol1'};
% %% Output:
% Singular mat file named seesion_name.mat

for i=1:numel(unitySessions)-1
    if i == 1
        [lickt,lickx,post,posx,trial,sp] = ...
            concatenateNPMatFiles(fullfile(data_dir,strcat(unitySessions{i},'.mat')),...
            fullfile(data_dir,strcat(unitySessions{i+1},'.mat')));

        save(fullfile(data_dir,strcat(session_name,'.mat')),'lickt','lickx','post','posx','trial','sp');
    else
        [lickt,lickx,post,posx,trial,sp] = ...
            concatenateNPMatFiles(fullfile(data_dir,strcat(session_name,'.mat')),...
            fullfile(data_dir,strcat(unitySessions{i+1},'.mat')));

        save(fullfile(data_dir,strcat(session_name,'.mat')),'lickt','lickx','post','posx','trial','sp');
    end
end

end




