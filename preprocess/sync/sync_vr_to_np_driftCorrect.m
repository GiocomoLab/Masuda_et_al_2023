function sync_vr_to_np_driftCorrect(data_dir, save_data_dir, session_name, tracklength, save_fig_bool)
fprintf('Start syncing session: %s @ %s\n',session_name,datestr(now,'mmmm dd, yyyy HH:MM:SS AM'));
 
addpath(genpath('C:\code\spikes'));
addpath(genpath('C:\code\npy-matlab'));
addpath(genpath('C:\code\AlexA_Library\'))
addpath(genpath('C:\code\MalcolmFxn\'))
addpath(genpath('C:\code\SiliconProbeCode\'))
 
% % specify tracklength
% tracklength = 400;
% 
% % location of data
% data_dir = 'Z:\giocomo\export\data\Projects\JohnKei_NPH3\G1\G1_190817_keicontrasttrack_ketamine1_g0';
% session_name = '190817_keicontrasttrack_baseline1';
% save_data_dir = 'Z:\giocomo\export\data\Projects\JohnKei_NPH3\data_giocomo\individualUnitySessionData';
img_save_data_dir = fullfile(save_data_dir,'syncImages');
 
[~,main_name]=fileparts(data_dir);
animal_name = strsplit(main_name,'_');
animal_name = animal_name{1};
NIDAQ_file = fullfile(data_dir,strcat(main_name,'_t0.nidq.bin'));
NIDAQ_config = fullfile(data_dir,strcat(main_name,'_t0.nidq.meta'));
spike_dir = fullfile(data_dir,strcat(main_name,'_imec0'));
 
%get the nidaq sample rate & get number of recorded nidaq channels
dat=textscan(fopen(NIDAQ_config),'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'niSampRate');
sync_sampling_rate=str2double(vals{loc});
 
loc2=contains(names,'nSavedChans');
n_channels_nidaq=str2double(vals{loc2});
 
% get neuropixels sync pulse times
fpNIDAQ=fopen(NIDAQ_file);
datNIDAQ=fread(fpNIDAQ,[n_channels_nidaq,Inf],'*int16');
fclose(fpNIDAQ);
syncDat=datNIDAQ(2,:)>1000;
 
%%
frame_times_np = find(abs(diff(syncDat))==1)+1;
frame_times_np = frame_times_np/sync_sampling_rate;
 
% read vr position data
formatSpec = '%f%f%f%f%f%[^\n\r]';
delimiter = '\t';
fid = fopen(fullfile(data_dir,strcat(session_name,'_position.txt')),'r');
dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fid);
vr_position_data = cat(2,dataArray{1:5});
%vr_position_data = vr_position_data(1:49334,:);
nu_entries = nnz(~isnan(vr_position_data(1,:)));
%vr_position_data=vr_position_data(49335:end,:);
vr_ttl=vr_position_data(:,nu_entries); %assuming TTL in last and timestamp in second last column
frame_times_vr=vr_position_data(:,nu_entries-1);
 
% read vr trial data
fid = fopen(fullfile(data_dir,strcat(session_name,'_trial_times.txt')),'r');
vr_trial_data = fscanf(fid, '%f', [4,inf])';
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
 
% set vr frame times to be the time of neuropixels pulses
% make sure the number of frames matches (can be off by one because of
% odd/even numbers of frames)
%%
tmp_diff=diff(frame_times_np);
[mm,step_idx]=find(tmp_diff>2); %2 
sess_length=diff([0 step_idx length(frame_times_np)]);
midpoint = ([0 step_idx] + [step_idx length(frame_times_np)])/2;
%step_idx=step_idx+1;
frametimes_nlOld = frame_times_np;
[~,ml]=min(abs(sess_length-numel(frame_times_vr)));
if length(mm)>=1
    
    if save_fig_bool
        figure(1);
        subplot(1,2,1)
        plot(tmp_diff)
        title(sprintf('found %d blocks',numel(step_idx+2)))
        hold on
        plot(step_idx,tmp_diff(step_idx),'ro')
        for im=1:numel(midpoint)
            text(midpoint(im),max(tmp_diff),sprintf('%d',sess_length(im)))
        end
    end
%     figure;
%     subplot(2,1,1)
%     plot(frame_times_np)
%     subplot(2,1,2)
%     plot(tmp_diff)
%     title(sprintf('found %d blocks',numel(step_idx+2)))
%     hold on
%     plot(step_idx,tmp_diff(step_idx),'ro')
%     
%     for im=1:numel(midpoint)
%     text(midpoint(im),max(tmp_diff),sprintf('%d',sess_length(im)))
%     end
    
%     sess=input(sprintf('Which session do you want to extract (suggesting %d)',ml));
    sess = ml;
    fprintf('Extracting session #%d \n',ml);
%     is_mismatch = input('Is this a MM or PB sesion [0/1]?');
    is_mismatch = 0;
 
    step_idx = [0 step_idx length(frame_times_np)];
    %frame_times_np=frame_times_np(ii+2:end);
    idx_start=step_idx(sess)+1;
    idx_stop = step_idx(sess+1);
    frame_times_np=frame_times_np(idx_start:idx_stop);
else
    is_mismatch=0;
end
%%
if abs(numel(frame_times_np) - numel(frame_times_vr)) > 1
    warning('number of sync pulses does not match number of frames.')
    fprintf('n NP: %d, n VR: %d \n',numel(frame_times_np),numel(frame_times_vr));
%     answer=input('Continue[0/1]?');
    answer = 1;
    if answer == 0
        error('Stopped')
    end
end
    idx=1:min(numel(frame_times_np),numel(frame_times_vr)); %use shorter index
    post = frame_times_np(idx)';
    vr_position_data=vr_position_data(idx,:);
    posx = vr_position_data(:,1);
    
    % transform lick times into neuropixels reference frame
    beta = [ones(size(post)) frame_times_vr(idx)]\post;
    lickt = beta(1) + lickt*beta(2);
%%
if save_fig_bool
    subplot(1,2,2)
    scatter(diff(post),diff(frame_times_vr(idx)),2,1:length(idx)-1)
    title(sprintf('Extracting session #%d \n',ml));
    saveas(figure(1), fullfile(img_save_data_dir,strcat(animal_name,'_',session_name,'.jpg')));
end
fprintf('Saved sync image @ %s\n',datestr(now,'mmmm dd, yyyy HH:MM:SS AM'));
 
%%
% load spike times
sp = loadKSdir(spike_dir);
 
%% CORRECT FOR DRIFT BETWEEN IMEC AND NIDAQ BOARDS
% TWO-PART CORRECTION
% 1. Get sync pulse times relative to NIDAQ and Imec boards.  
% 2. Quantify difference between the two sync pulse times and correct in
% spike.st. 
% PART 1: GET SYNC TIMES RELATIVE TO EACH BOARD
% We already loaded most of the NIDAQ data above. Here, we access the sync
% pulses used to sync Imec and NIDAQ boards together. The times a pulse is
% emitted and registered by the NIDAQ board are stored in syncDatNIDAQ below.
syncDatNIDAQ=datNIDAQ(1,:)>1000;
% convert NIDAQ sync data into time data by dividing by the sampling rate
ts_NIDAQ = strfind(syncDatNIDAQ,[0 1])/sync_sampling_rate; 
% ts_NIDAQ: these are the sync pulse times relative to the NIDAQ board
% Now, we do the same, but from the perspective of the Imec board. 
LFP_config = dir(fullfile(spike_dir,'*.lf.meta'));
LFP_config = fullfile(LFP_config.folder,LFP_config.name);
LFP_file = dir(fullfile(spike_dir,'*.lf.bin'));
LFP_file = fullfile(LFP_file.folder,LFP_file.name);
dat=textscan(fopen(LFP_config),'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'imSampRate');
lfp_sampling_rate=str2double(vals{loc});
% for loading only a portion of the LFP data
fpLFP = fopen(LFP_file);
fseek(fpLFP, 0, 'eof'); % go to end of file
fpLFP_size = ftell(fpLFP); % report size of file
fpLFP_size = fpLFP_size/(2*384); 
fclose(fpLFP);
% get the sync pulse times relative to the Imec board
fpLFP=fopen(LFP_file);
fseek(fpLFP,384*2,0);
ftell(fpLFP);
datLFP=fread(fpLFP,[1,round(fpLFP_size/4)],'*int16',384*2); % this step used to take forever
fclose(fpLFP);
syncDatLFP=datLFP(1,:)>10; 
ts_LFP = strfind(syncDatLFP,[0 1])/lfp_sampling_rate;
% ts_LFP: these are the sync pulse times relative to the Imec board
% PART 2: TIME CORRECTION
lfpNIDAQdif = ts_LFP - ts_NIDAQ(1:size(ts_LFP, 2)); % calculate the difference between the sync pulse times
fit = polyfit(ts_LFP, lfpNIDAQdif, 1); % linear fit 
correction_slope = fit(1); % this is the amount of drift we get per pulse (that is, per second)
% save the old, uncorrected data as sp.st_uncorrected and save the new,
% corrected data as sp.st (as many of your analyses are using sp.st).
sp.st_uncorrected = sp.st; % save uncorrected spike times (st)
st_corrected = sp.st - sp.st * correction_slope; % in two steps to avoid confusion
sp.st = st_corrected; % overwrite the old sp.st
%%
% shift everything to start at zero
offset = post(1);
post = post - offset;
lickt = lickt - offset;
sp.st = sp.st - offset;
sp.vr_session_offset = offset;
 
% resample position to have constant time bins
posx = interp1(post,posx,(0:0.02:max(post))');
vr_data_downsampled=interp1(post,vr_position_data,(0:0.02:max(post)));
post = (0:0.02:max(post))';
posx([false;diff(posx)<-2])=round(posx([false;diff(posx)<-2])/tracklength)*tracklength; % handle teleports
 
% compute trial number for each time bin
trial = [1; cumsum(diff(posx)<-100)+1]; 
 
% throw out bins after the last trial
if is_mismatch
    num_trials = max(trial);
end
keep = trial<=num_trials;
trial = trial(keep);
posx = posx(keep);
post = post(keep);
 
% throw out licks before and after session
keep = lickt>min(post) & lickt<max(post);
lickx = lickx(keep);
lickt = lickt(keep);
 
% cut off all spikes before and after vr session
keep = sp.st >= 0 & sp.st <= post(end);
sp.st = sp.st(keep);
sp.spikeTemplates = sp.spikeTemplates(keep);
sp.clu = sp.clu(keep);
sp.tempScalingAmps = sp.tempScalingAmps(keep);
 
% save processed data
if is_mismatch == 0
save(fullfile(save_data_dir,strcat(animal_name,'_',session_name,'.mat')),'sp','post','posx','lickt','lickx','trial','trial_contrast','trial_gain');
%save(fullfile(data_dir,strcat(session_name,'.mat')),'sp','post','posx','lickt','lickx','trial','trial_contrast','trial_gain');
 
else
 true_speed= vr_data_downsampled(:,2);
 mismatch_trigger = vr_data_downsampled(:,3);
 save(fullfile(save_data_dir,strcat(animal_name,'_',session_name,'.mat')),'sp','post','posx','lickt','lickx','trial','trial_contrast','trial_gain','true_speed','mismatch_trigger');
end
 
fprintf('Finished syncing session: %s @ %s\n',session_name,datestr(now,'mmmm dd, yyyy HH:MM:SS AM'));

