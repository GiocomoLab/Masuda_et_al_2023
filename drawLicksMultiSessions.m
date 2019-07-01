function drawLicksMultiSessions(data_dir)

%% Multiple Session Stitch
close all;
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
% where to find data and save images
% data_dir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190619_keicontrasttrack_ketamine1_g0';
[~,name,~] = fileparts(data_dir);
session_root = strsplit(name,'track');
session_name = {strcat(session_root{1},'track_baseline+cntrlinjx+ketamine')};
session_name1  = {strcat(session_root{1},'track_baseline1')};
session_name2 = {strcat(session_root{1},'track_controlinjx1')};
session_name3 = {strcat(session_root{1},'track_ketamine1')};

% load data
fprintf('session %d/%d: %s\n','1','1',session_name{1});
load(fullfile(data_dir,strcat(session_name3{1},'.mat')),'lickt','lickx','post','posx','sp','trial');%ketamine 2
lickt2 = lickt;lickx2 = lickx;post2 = post;posx2=posx;sp2 = sp;trial2 = trial;
load(fullfile(data_dir,strcat(session_name2{1},'.mat')),'lickt','lickx','post','posx','sp','trial');%postsaline 1
lickt1 = lickt;lickx1 = lickx;post1 = post;posx1=posx;sp1 = sp;trial1 = trial;
load(fullfile(data_dir,strcat(session_name1{1},'.mat')),'lickt','lickx','post','posx','sp','trial'); %presaline 0

[~,~,lick_idx] = histcounts(lickt,post);
[~,~,lick_idx1] = histcounts(lickt1,post1);
[~,~,lick_idx2] = histcounts(lickt2,post2);
%%
close all;
clear g;
h = figure('Position',[100 100 1200 600]);
g(1,1) = gramm('y',vertcat(trial(lick_idx)-100,trial1(lick_idx1)-50,trial2(lick_idx2)),'x',vertcat(posx(lick_idx),posx1(lick_idx1),posx2(lick_idx2)));
g(1,1).geom_point();
g(1,1).set_title(session_name);
g(1,1).set_names('x','Lick Position in VR','y','trial number');

g(1,2) = gramm('y',vertcat(trial(lick_idx)-100,trial1(lick_idx1)-50,trial2(lick_idx2)),'x',vertcat(posx(lick_idx),posx1(lick_idx1),posx2(lick_idx2)));
g(1,2).stat_bin();
g(1,2).set_title('Binned Licks by Position');
g(1,2).set_names('x','Lick Position in VR','y','Lick Number');

g(1,3) = gramm('x',vertcat(trial(lick_idx)-100,trial1(lick_idx1)-50,trial2(lick_idx2)),'y',vertcat(posx(lick_idx),posx1(lick_idx1),posx2(lick_idx2)));
g(1,3).stat_bin('nbins',30,'geom','line');
g(1,3).set_title('Binned Licks by Trial');
g(1,3).set_names('x','Trial Number','y','Number of Licks per 10 Trials');
g.draw();

%%
% make image dir if it doesn't exist
image_save_dir = fullfile('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/images/behavior_plots');
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
saveas(h,fullfile(image_save_dir,session_name{1}),'png');
fprintf(image_save_dir);

end