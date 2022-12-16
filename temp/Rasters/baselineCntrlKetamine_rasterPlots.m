% plots rasters for all cells in 1 sessions

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
data_dir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCN1/HCN1_190623_keicontrasttrack_ketamine1_g0';
session_name = {'HCN1_190623_keicontrasttrack_baseline+cntrlinjx+ketamine'};
session_name1  = {'HCN1_190623_keicontrasttrack_baseline1'};
session_name2 = {'HCN1_190623_keicontrasttrack_controlinjx1'};
session_name3 = {'HCN1_190623_keicontrasttrack_ketamine2'};
% load data
fprintf('session %d/%d: %s\n','1','1',session_name{1});
load(fullfile(data_dir,strcat(session_name3{1},'.mat')),'lickt','lickx','post','posx','sp','trial');%ketamine 2
lickt2 = lickt;lickx2 = lickx;post2 = post;posx2=posx;sp2 = sp;trial2 = trial;
load(fullfile(data_dir,strcat(session_name2{1},'.mat')),'lickt','lickx','post','posx','sp','trial');%postsaline 1
lickt1 = lickt;lickx1 = lickx;post1 = post;posx1=posx;sp1 = sp;trial1 = trial;
load(fullfile(data_dir,strcat(session_name1{1},'.mat')),'lickt','lickx','post','posx','sp','trial'); %presaline 0


cells_to_plot = sp.cids(sp.cgs==2); % all good cells
cells_to_plot1 = sp1.cids(sp1.cgs==2); % all good cells
cells_to_plot2 = sp2.cids(sp2.cgs==2); % all good cells

% make image dir if it doesn't exist
image_save_dir = strcat('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/images/',...
    session_name{1},'/pretty_rasters/');
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
%%
% compute some useful information (like spike depths)
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

[spikeAmps1, spikeDepths1, templateYpos1, tempAmps1, tempsUnW1, tempDu1r, tempPeakWF1] = ...
    templatePositionsAmplitudes(sp1.temps, sp1.winv, sp1.ycoords, sp1.spikeTemplates, sp1.tempScalingAmps);

[spikeAmps2, spikeDepths2, templateYpos2, tempAmps2, tempsUnW2, tempDu2r, tempPeakWF2] = ...
    templatePositionsAmplitudes(sp2.temps, sp2.winv, sp2.ycoords, sp2.spikeTemplates, sp2.tempScalingAmps);

% get spike depths
spike_depth = nan(numel(cells_to_plot),1);
for k = 1:numel(cells_to_plot)
    spike_depth(k) = median(spikeDepths(sp.clu==cells_to_plot(k)));
end

% get spike depths 1 
spike_depth1 = nan(numel(cells_to_plot1),1);
for k = 1:numel(cells_to_plot1)
    spike_depth1(k) = median(spikeDepths1(sp1.clu==cells_to_plot1(k)));
end

% get spike depths 2
spike_depth2 = nan(numel(cells_to_plot2),2);
for k = 2:numel(cells_to_plot2)
    spike_depth2(k) = median(spikeDepths2(sp2.clu==cells_to_plot2(k)));
end
% sort cells_to_plot by spike_depth (descending)
[spike_depth,sort_idx] = sort(spike_depth,'descend');
cells_to_plot = cells_to_plot(sort_idx);

% sort cells_to_plot by spike_depth (descending) 1 
[spike_depth1,sort_idx1] = sort(spike_depth1,'descend');
cells_to_plot1 = cells_to_plot1(sort_idx1);

% sort cells_to_plot by spike_depth (descending)
[spike_depth2,sort_idx2] = sort(spike_depth2,'descend');
cells_to_plot2 = cells_to_plot2(sort_idx2);
%% make raster plots for all cells
h = figure('Position',[100 100 160 500]); hold on;

for k = 1:numel(cells_to_plot)
    cla;
    
    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

    % get spike times and index into post
    spike_t = sp.st(sp.clu==cells_to_plot(k));
    [~,~,spike_idx] = histcounts(spike_t,post);
    
    % get spike times and index into post1
    spike_t1 = sp1.st(sp1.clu==cells_to_plot(k));
    [~,~,spike_idx1] = histcounts(spike_t1,post1);
    
    % get spike times and index into post2
    spike_t2 = sp2.st(sp2.clu==cells_to_plot(k));
    [~,~,spike_idx2] = histcounts(spike_t2,post2);
    
    % [~,~,spike_idx1_20min] = histcounts(spike_t1(spike_t1<1200),post1); %color first 20min after ketamine red
    
    % plot spike raster
    presalineTrialNum = 50;
    postsalineTrialNum = 50;
    
    plot(posx(spike_idx),trial(spike_idx)-(presalineTrialNum+postsalineTrialNum),'k.');
    plot(posx1(spike_idx1),trial1(spike_idx1)-postsalineTrialNum,'k.'); %b.
    plot(posx2(spike_idx2),trial2(spike_idx2),'k.'); %r.
    % plot(posx1(spike_idx1_20min),trial1(spike_idx1_20min),'r.');
    xlim([params.TrackStart params.TrackEnd]);
    ylim([-(max(trial)+max(trial1)) max(trial2)+1]);
    %ylim([-max(trial+1) max(trial1)+1]);
    title(sprintf('c%d, d=%d',cells_to_plot(k),round(spike_depth(k))));
    % xticks(''); yticks('');
    

    % save fig
    if save_figs
        saveas(h,fullfile(image_save_dir,sprintf('%d.png',k)),'png');
        saveas(h,fullfile(image_save_dir,sprintf('%d.pdf',k)),'pdf');
    end
    

end