function drawTrainingLicks(data_dir,session_name)

% location of data
% data_dir = 'F:\H3\npH3_0401_gain_1_g0';
[~,main_name]=fileparts(data_dir);
% session_name = '0401_dark_3';

% read vr position data
formatSpec = '%f%f%f%f%f%[^\n\r]';
delimiter = '\t';
fid = fopen(fullfile(data_dir,strcat(session_name,'_position.txt')),'r');
dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fid);
vr_position_data = cat(2,dataArray{1:5});
nu_entries = nnz(~isnan(vr_position_data(1,:)));
vr_ttl=vr_position_data(:,nu_entries); %assuming TTL in last and timestamp in second last column
frame_times_vr=vr_position_data(:,nu_entries-1);

posx = vr_position_data(:,nu_entries-2);
post = frame_times_vr;
% compute trial number for each time bin
trial = [1; cumsum(diff(posx)<-100)+1];


% read vr trial data
fid = fopen(fullfile(data_dir,strcat(session_name,'_trial_times.txt')),'r');
vr_trial_data =  fscanf(fid, '%f', [4,inf])';
fclose(fid);
trial_contrast = [100; vr_trial_data(:,2)];
trial_gain = [1; vr_trial_data(:,3)];
num_trials = numel(trial_gain);

% read vr licking data
fid = fopen(fullfile(data_dir,strcat(session_name,'_licks.txt')),'r');
vr_lick_data = fscanf(fid, '%f', [2,inf])';
fclose(fid);
lickx = vr_lick_data(:,1);
lickt = vr_lick_data(:,2);

%%

[~,~,lick_idx] = histcounts(lickt,post);

lickAccuracyByTrial = zeros(1,max(max(trial)));
for i = 1:max(max(trial))
    trialLicks = lickx(trial(lick_idx) == i);
    goodLicks = sum(trialLicks<25) + sum(trialLicks>max(posx)-25); 
    if trialLicks ~= 0
        lickAccuracyByTrial(i) = goodLicks/numel(trialLicks);
    else
        lickAccuracyByTrial(i) = 0.0;
    end
end


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

g(1,4) = gramm('x',1:max(max(trial)),'y',lickAccuracyByTrial);
g(1,4).stat_smooth();
g(1,4).set_title('Lick Accuracy - Smoothed');
g(1,4).set_names('x','Trial Number','y','Lick Accuracy (% within 50cm of target)');
g.draw();



end