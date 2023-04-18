function [start_time, end_time] = get_vr_start_end(data_dir, session_name)

[~,main_name]=fileparts(data_dir);
NIDAQ_file = fullfile(data_dir,strcat(main_name,'_t0.nidq.bin'));
NIDAQ_config = fullfile(data_dir,strcat(main_name,'_t0.nidq.meta'));

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
nu_entries = nnz(~isnan(vr_position_data(1,:)));
frame_times_vr=vr_position_data(:,nu_entries-1);

% set vr frame times to be the time of neuropixels pulses
% make sure the number of frames matches (can be off by one because of
% odd/even numbers of frames)
tmp_diff=diff(frame_times_np);
[mm,step_idx]=find(tmp_diff>2); %2
sess_length=diff([0 step_idx length(frame_times_np)]);
[~,ml]=min(abs(sess_length-numel(frame_times_vr)));
if length(mm)>=1
    sess = ml;
    step_idx = [0 step_idx length(frame_times_np)];
    idx_start=step_idx(sess)+1;
    idx_stop = step_idx(sess+1);
    frame_times_np=frame_times_np(idx_start:idx_stop);
end

idx=1:min(numel(frame_times_np),numel(frame_times_vr)); %use shorter index
start_time = frame_times_np(1);
end_time = frame_times_np(idx(end));