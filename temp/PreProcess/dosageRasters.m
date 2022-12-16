%% params
% make sure paths are correct
restoredefaultpath
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
% some params
params = readtable('UniversalParams.xlsx');
save_figs = true;
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;

% where to find data and save images
data_dir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/A3/2019-05-20_09-46-57';
session_name = {'B3_0520_contrasttrack_dosage'};
session_name0  = {'0520_contrasttrack_baseline1'};
session_name1  = {'0520_contrasttrack_5mgkg_ketamine1'};
session_name2 = {'0520_contrasttrack_10mgkg_ketamine1'};
session_name3 = {'0520_contrasttrack_25mgkg_ketamine1'};
session_name4 = {'0520_contrasttrack_50mgkg_ketamine1'};
session_name5 = {'0520_contrasttrack_100mgkg_ketamine1'};

% load data
fprintf('session %d/%d: %s\n','1','1',session_name{1});
load(fullfile(data_dir,strcat(session_name5{1},'.mat')),'lickt','lickx','post','posx','sp','trial');%ketamine 2
lickt5 = lickt;lickx5 = lickx;post5 = post;posx5=posx;sp5 = sp;trial5 = trial;
load(fullfile(data_dir,strcat(session_name4{1},'.mat')),'lickt','lickx','post','posx','sp','trial');%ketamine 2
lickt4 = lickt;lickx4 = lickx;post4 = post;posx4=posx;sp4 = sp;trial4 = trial;
load(fullfile(data_dir,strcat(session_name3{1},'.mat')),'lickt','lickx','post','posx','sp','trial');%ketamine 2
lickt3 = lickt;lickx3 = lickx;post3 = post;posx3=posx;sp3 = sp;trial3 = trial;
load(fullfile(data_dir,strcat(session_name2{1},'.mat')),'lickt','lickx','post','posx','sp','trial');%ketamine 2
lickt2 = lickt;lickx2 = lickx;post2 = post;posx2=posx;sp2 = sp;trial2 = trial;
load(fullfile(data_dir,strcat(session_name1{1},'.mat')),'lickt','lickx','post','posx','sp','trial');%postsaline 1
lickt1 = lickt;lickx1 = lickx;post1 = post;posx1=posx;sp1 = sp;trial1 = trial;
load(fullfile(data_dir,strcat(session_name0{1},'.mat')),'lickt','lickx','post','posx','sp','trial'); %presaline 0
lickt0 = lickt;lickx0 = lickx;post0 = post;posx0=posx;sp0 = sp;trial0 = trial;

cells_to_plot0 = sp0.cids(sp0.cgs==2); % all good cells
cells_to_plot1 = sp1.cids(sp1.cgs==2); % all good cells
cells_to_plot2 = sp2.cids(sp2.cgs==2); % all good cells
cells_to_plot3 = sp3.cids(sp3.cgs==2); % all good cells
cells_to_plot4 = sp4.cids(sp4.cgs==2); % all good cells
cells_to_plot5 = sp5.cids(sp5.cgs==2); % all good cells

% make image dir if it doesn't exist
image_save_dir = strcat('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/images/',...
    session_name{1},'/pretty_rasters/');
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
%%
% compute some useful information (like spike depths)
[spikeAmps0, spikeDepths0, templateYpos0, tempAmps0, tempsUnW0, tempDur0, tempPeakWF0] = ...
    templatePositionsAmplitudes(sp0.temps, sp0.winv, sp0.ycoords, sp0.spikeTemplates, sp0.tempScalingAmps);

[spikeAmps1, spikeDepths1, templateYpos1, tempAmps1, tempsUnW1, tempDur1, tempPeakWF1] = ...
    templatePositionsAmplitudes(sp1.temps, sp1.winv, sp1.ycoords, sp1.spikeTemplates, sp1.tempScalingAmps);

[spikeAmps2, spikeDepths2, templateYpos2, tempAmps2, tempsUnW2, tempDur2, tempPeakWF2] = ...
    templatePositionsAmplitudes(sp2.temps, sp2.winv, sp2.ycoords, sp2.spikeTemplates, sp2.tempScalingAmps);

[spikeAmps3, spikeDepths3, templateYpos3, tempAmps3, tempsUnW3, tempDur3, tempPeakWF3] = ...
    templatePositionsAmplitudes(sp3.temps, sp3.winv, sp3.ycoords, sp3.spikeTemplates, sp3.tempScalingAmps);

[spikeAmps4, spikeDepths4, templateYpos4, tempAmps4, tempsUnW4, tempDur4, tempPeakWF4] = ...
    templatePositionsAmplitudes(sp4.temps, sp4.winv, sp4.ycoords, sp4.spikeTemplates, sp4.tempScalingAmps); 

[spikeAmps5, spikeDepths5, templateYpos5, tempAmps5, tempsUnW5, tempDur5, tempPeakWF5] = ...
    templatePositionsAmplitudes(sp5.temps, sp5.winv, sp5.ycoords, sp5.spikeTemplates, sp5.tempScalingAmps); 


% get spike depths 0
spike_depth0 = nan(numel(cells_to_plot0),1);
for k = 1:numel(cells_to_plot0)
    spike_depth0(k) = median(spikeDepths0(sp0.clu==cells_to_plot0(k)));
end

% get spike depths 1 
spike_depth1 = nan(numel(cells_to_plot1),1);
for k = 1:numel(cells_to_plot1)
    spike_depth1(k) = median(spikeDepths1(sp1.clu==cells_to_plot1(k)));
end

% get spike depths 2
spike_depth2 = nan(numel(cells_to_plot2),1);
for k = 2:numel(cells_to_plot2)
    spike_depth2(k) = median(spikeDepths2(sp2.clu==cells_to_plot2(k)));
end

% get spike depths 3
spike_depth3 = nan(numel(cells_to_plot3),1);
for k = 2:numel(cells_to_plot3)
    spike_depth3(k) = median(spikeDepths3(sp3.clu==cells_to_plot3(k)));
end

% get spike depths 4
spike_depth4 = nan(numel(cells_to_plot4),1);
for k = 2:numel(cells_to_plot4)
    spike_depth4(k) = median(spikeDepths4(sp4.clu==cells_to_plot4(k)));
end

% get spike depths 5
spike_depth5 = nan(numel(cells_to_plot5),1);
for k = 2:numel(cells_to_plot5)
    spike_depth5(k) = median(spikeDepths5(sp5.clu==cells_to_plot4(k)));
end


% sort cells_to_plot by spike_depth (descending) 0
[spike_depth0,sort_idx0] = sort(spike_depth0,'descend');
cells_to_plot0 = cells_to_plot0(sort_idx0);

% sort cells_to_plot by spike_depth (descending) 1 
[spike_depth1,sort_idx1] = sort(spike_depth1,'descend');
cells_to_plot1 = cells_to_plot1(sort_idx1);

% sort cells_to_plot by spike_depth (descending)
[spike_depth2,sort_idx2] = sort(spike_depth2,'descend');
cells_to_plot2 = cells_to_plot2(sort_idx2);

% sort cells_to_plot by spike_depth (descending)
[spike_depth3,sort_idx3] = sort(spike_depth3,'descend');
cells_to_plot3 = cells_to_plot3(sort_idx3);

% sort cells_to_plot by spike_depth (descending)
[spike_depth4,sort_idx4] = sort(spike_depth4,'descend');
cells_to_plot4 = cells_to_plot4(sort_idx4);

% sort cells_to_plot by spike_depth (descending)
[spike_depth5,sort_idx5] = sort(spike_depth5,'descend');
cells_to_plot5 = cells_to_plot5(sort_idx5);
%% make raster plots for all cells
h = figure('Position',[100 100 160 500]); hold on;

for k = 1:numel(cells_to_plot)
    cla;
    
    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

    % get spike times and index into post1
    spike_t0 = sp0.st(sp0.clu==cells_to_plot(k));
    [~,~,spike_idx0] = histcounts(spike_t0,post0);

    % get spike times and index into post1
    spike_t1 = sp1.st(sp1.clu==cells_to_plot(k));
    [~,~,spike_idx1] = histcounts(spike_t1,post1);
    
    % get spike times and index into post2
    spike_t2 = sp2.st(sp2.clu==cells_to_plot(k));
    [~,~,spike_idx2] = histcounts(spike_t2,post2);

        % get spike times and index into post2
    spike_t3 = sp3.st(sp3.clu==cells_to_plot(k));
    [~,~,spike_idx3] = histcounts(spike_t3,post3);

    % get spike times and index into post2
    spike_t4 = sp4.st(sp4.clu==cells_to_plot(k));
    [~,~,spike_idx4] = histcounts(spike_t4,post4);
    
    % get spike times and index into post2
    spike_t5 = sp5.st(sp5.clu==cells_to_plot(k));
    [~,~,spike_idx5] = histcounts(spike_t5,post5);
    % [~,~,spike_idx1_20min] = histcounts(spike_t1(spike_t1<1200),post1); %color first 20min after ketamine red
    
    % plot spike raster
    baseline = 50;
    postsalineTrialNum = 100;
    
    plot(posx0(spike_idx0),trial0(spike_idx0),'k.'); %r.
    plot(posx1(spike_idx1),trial1(spike_idx1)+50,'k.'); %b.
    plot(posx2(spike_idx2),trial2(spike_idx2)+150,'k.'); %r.
    plot(posx3(spike_idx3),trial3(spike_idx3)+250,'k.'); %b.
    plot(posx4(spike_idx4),trial4(spike_idx4)+350,'k.'); %r.
    plot(posx5(spike_idx5),trial5(spike_idx5)+450,'k.'); %r.
    % plot(posx1(spike_idx1_20min),trial1(spike_idx1_20min),'r.');
    xlim([params.TrackStart params.TrackEnd]);
    %ylim([-(max(trial)+max(trial1)) max(trial2)+1]);
    ylim([-1 551]);
    title(sprintf('c%d, d=%d',cells_to_plot(k),round(spike_depth(k))));
    % xticks(''); yticks('');
    

    % save fig
    if save_figs
        saveas(h,fullfile(image_save_dir,sprintf('%d.png',k)),'png');
        saveas(h,fullfile(image_save_dir,sprintf('%d.pdf',k)),'pdf');
    end
    

end