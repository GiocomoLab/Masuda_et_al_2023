function lickAccuracyByTrial = drawLicksSingleSessions(data_dir,session_name)
% Inputs: data directory and session name
% Assumes reward is at end of track. 'Accurate' lick = within 50cm of reward zone
% Image with lick raster plot + associated quantification plots

addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
% where to find data and save images
% data_dir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190619_keicontrasttrack_ketamine1_g0';
[~,name,~] = fileparts(data_dir);
% session_root = strsplit(name,'track');
%session_name = {strcat(session_root{1},'track_baseline+cntrlinjx+ketamine')};


% load data
fprintf('session: %s\n',session_name);
load(fullfile(data_dir,strcat(session_name,'.mat')),'lickt','lickx','post','posx','sp','trial'); %presaline 0

[~,~,lick_idx] = histcounts(lickt,post);

lickAccuracyByTrial = zeros(1,max(trial));
for i = 1:max(trial)
    trialLicks = lickx(trial(lick_idx) == i);
    goodLicks = sum(trialLicks<25) + sum(trialLicks>max(posx)-25); 
    if trialLicks ~= 0
        lickAccuracyByTrial(i) = goodLicks/numel(trialLicks);
    else
        lickAccuracyByTrial(i) = 0.0;
    end
end


%%
close all;
clear g;
h = figure('Position',[100 100 1800 600]);
g(1,1) = gramm('y',trial(lick_idx),'x',posx(lick_idx));
g(1,1).geom_point();
g(1,1).set_title(session_name);
g(1,1).set_names('x','Lick Position in VR','y','trial number');

g(1,2) = gramm('y',trial(lick_idx),'x',posx(lick_idx));
g(1,2).stat_bin();
g(1,2).set_title('Binned Licks by Position');
g(1,2).set_names('x','Lick Position in VR','y','Lick Number');

g(1,3) = gramm('x',trial(lick_idx),'y',posx(lick_idx));
g(1,3).stat_bin('nbins',30,'geom','line');
g(1,3).set_title('Binned Licks by Trial');
g(1,3).set_names('x','Trial Number','y','Number of Licks per 10 Trials');

g(1,4) = gramm('x',1:max(trial),'y',lickAccuracyByTrial);
g(1,4).geom_line();
g(1,4).set_title('Lick Accuracy');
g(1,4).set_names('x','Trial Number','y','Lick Accuracy (% within 50cm of target)');
g.draw();

%%
% make image dir if it doesn't exist
image_save_dir = fullfile('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/images/behavior_plots');
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
saveas(h,fullfile(image_save_dir,session_name),'png');
fprintf(strcat(image_save_dir,'\n'));


end