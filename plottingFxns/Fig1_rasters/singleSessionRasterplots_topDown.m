function singleSessionRasterplots_topDown(data_dir,session_name, trackLength, plotWidth, plotHeight)
% plots rasters for all cells in 1 session given data directory and name of
% mat file that has neuropixel data synched with unity information
% Requires Access ot the serer
% Inputs:
%   data_dir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190625_keicontrasttrack_propofol1_g0';
%   session_name = 'G4_190625_keicontrasttrack_baseline+cntrlinjx+propofol';
%   trackLength = 400;
%   plotWidth = 160
%   plotHeight = 500
%   drugFlag = 100 (0 for regular plotting; 100 for plotting the drug injx = 0
% Output: 
%   Rasterplots saved to Oak


    animalName = extractBefore(session_name,'_');
    sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
    
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


    % load data
    fprintf('session: %s\n',session_name);
    load(fullfile(data_dir,strcat(session_name,'.mat')),'lickt','lickx','post','posx','sp','trial'); %presaline 0


    cells_to_plot = sp.cids(sp.cgs==2); % all good cells


    % make image dir if it doesn't exist
    image_save_dir = strcat('/Users/KeiMasuda/Desktop/fkm_analysis/topDownRasterplots/');
    if exist(image_save_dir,'dir')~=7
        mkdir(image_save_dir);
    end
    
    % compute some useful information (like spike depths)
    [~, spikeDepths, ~, ~, ~, ~, ~] = ...
        templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);


    % get spike depths
    spike_depth = nan(numel(cells_to_plot),1);
    for k = 1:numel(cells_to_plot)
        spike_depth(k) = median(spikeDepths(sp.clu==cells_to_plot(k)));
    end
    
    % Calculate 
    % TODO

    % sort cells_to_plot by spike_depth (descending)
    [spike_depth,sort_idx] = sort(spike_depth,'descend');
    cells_to_plot = cells_to_plot(sort_idx);

    % make raster plots for all cells
    h = figure('Position',[100 100 plotWidth plotHeight]); hold on; % changed from 160 to 320 for repeating, also from 500 to 166 for 50 trials
    axis ij
    
    for k = 1:numel(cells_to_plot)
        cla;

        fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));

        % get spike times and index into post
        spike_t = sp.st(sp.clu==cells_to_plot(k));
        [~,~,spike_idx] = histcounts(spike_t,post);

        plot(posx(spike_idx),trial(spike_idx),'k.');

        ylim([0 max(trial)+1]);

        xlim([params.TrackStart trackLength]);
%         title(sprintf('c%d, d=%d',cells_to_plot(k),round(spike_depth(k))));
        title(sprintf('%s-%s,\ncell:%d',animalName,sessionDate,cells_to_plot(k)));
        % xticks(''); yticks('');


        % save fig
        if save_figs
            saveas(h,fullfile(image_save_dir,sprintf('%s%s%s%s%d.png',animalName,'_',sessionDate,'_',k)),'png');
        end


    end

end