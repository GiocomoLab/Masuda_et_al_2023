% plots rasters for all cells in 1 sessions
% run first on normal trials, comment out spike_depth part and sorting
% based on spike_depth and cell_to_plot , then run on meth trials

% params
% make sure paths are correct
restoredefaultpath
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));

addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));

% some params
params = readtable('UniversalParams.xlsx');
save_figs = true;
xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;

% where to find data and save images
data_dir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/F3/F3_190625_johnrepeatingtrack_train1_g0';
session_name = {'F3_190625_johnrepeatingtrack_meth2'}; % output
session_name1  = {'F3_190625_johnrepeatingtrack_meth2'}; % looking for .mat data file (input)
trackLength = 640;



% load data
fprintf('session %d/%d: %s\n','1','1',session_name{1});
load(fullfile(data_dir,strcat(session_name1{1},'.mat')),'lickt','lickx','post','posx','sp','trial'); %presaline 0


% cells_to_plot = sp.cids(sp.cgs==2); % all good cells (COMMENT THIS OUT)


% make image dir if it doesn't exist
image_save_dir = strcat('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/images/',...
    session_name{1},'/pretty_rasters/');
if exist(image_save_dir,'dir')~=7
    mkdir(image_save_dir);
end
%
% compute some useful information (like spike depths)
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);

% COMMENT THIS OUT
% % get spike depths
% spike_depth = nan(numel(cells_to_plot),1);
% for k = 1:numel(cells_to_plot)
%     spike_depth(k) = median(spikeDepths(sp.clu==cells_to_plot(k)));
% end
% 
% % sort cells_to_plot by spike_depth (descending)
% [spike_depth,sort_idx] = sort(spike_depth,'descend');
% cells_to_plot = cells_to_plot(sort_idx);

% make raster plots for all cells
h = figure('Position',[100 100 320 80]); hold on; % changed from 160 to 320 for repeating, also from 500 to 166 for 50 trials

for k = 1:numel(cells_to_plot)
    cla;
    
    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

    % get spike times and index into post
    spike_t = sp.st(sp.clu==cells_to_plot(k));
    [~,~,spike_idx] = histcounts(spike_t,post);

    
    plot(posx(spike_idx),trial(spike_idx),'k.');
    
    xlim([params.TrackStart trackLength]);
    ylim([0 max(trial)+1]);
    title(sprintf('c%d, d=%d',cells_to_plot(k),round(spike_depth(k))));
    % xticks(''); yticks('');
    

    % save fig
    if save_figs
        saveas(h,fullfile(image_save_dir,sprintf('%d.png',k)),'png');
        saveas(h,fullfile(image_save_dir,sprintf('%d.pdf',k)),'pdf');
    end
    
    
end