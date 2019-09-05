        
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));

addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));
% some params
params = readtable('UniversalParams.xlsx');
p = params;
save_figs = true;
%%

sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/*.mat');
mtrx_sesh = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices/*.mat');
%% Remove strange sessions
remove = {'AA','B1','B3','E1','E2','F3'}; %check for this in session name and remove
for z= 1:numel(remove)
    idx = ~cellfun('isempty',strfind({sessions.name},remove{z}));
    sessions(idx) = [];
end

%%
totalCell = 1;
allRho = [];
for n = 1:numel(sessions)
    try
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        session_name = sessions(n).name(1:end-4);
        
        animalName = extractBefore(session_name,'_');
        sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
        
        fprintf(strcat('\nProcessing:', session_name,'\n'));
        trackLength = 400;
        trackEnd = trackLength;
        plotWidth = 160;
        plotHeight = 500;
        preDrugTrials = 100;
        trials_corrTemplate = 50;
        
        xbincent = params.TrackStart+params.SpatialBin/2:params.SpatialBin:params.TrackEnd-params.SpatialBin/2;


        % load data
        fprintf('session: %s\n',session_name);
        load(matPath,'lickt','lickx','post','posx','sp','trial'); 


        cells_to_plot = sp.cids(sp.cgs==2); % all good cells

%%
        % make image dir if it doesn't exist
        image_save_dir = strcat('/Users/KeiMasuda/Desktop/fkm_analysis/img/pretty_rasters/');
        if exist(image_save_dir,'dir')~=7
            mkdir(image_save_dir);
        end

        % compute some useful information (like spike depths)
        [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
            templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);


        % get spike depths
        spike_depth = nan(numel(cells_to_plot),1);
        for k = 1:numel(cells_to_plot)
            spike_depth(k) = median(spikeDepths(sp.clu==cells_to_plot(k)));
        end


        % sort cells_to_plot by spike_depth (descending)
        [spike_depth,sort_idx] = sort(spike_depth,'descend');
        cells_to_plot = cells_to_plot(sort_idx);
        
        %% Calculate speed and filter out stationary periods (speed<2cm)
        speed = calcSpeed(posx, p);
        [trial,posx,post] = speedFilterData(trial,posx, post, speed);
        %% Preallocate and store data for all cell's firing rates
        nCells = size(cells_to_plot, 2);
        spatialBins = trackLength/params.SpatialBin; 
        all_fr = nan(nCells, max(trial),spatialBins); % preallocate matrix of cells' firing rates across spatial bins

        all_corrmatrix = nan(nCells, max(trial), max(trial)); % preallocate matrix of cells' correlations between trials
        % all_corrblock = nan(nCells, floor(max(trial)/trials_per_block), floor(max(trial)/trials_per_block)); % preallocate matrix of cells' correlations every 50 trials
%         all_waveforms= nan(nCells, size(waveforms,2)); 
        all_cellCorrScore = nan(nCells, numel(trials_corrTemplate:max(trial)));
        all_correlationScore = nan(nCells, 2);
        fprintf('Calculating firing rate for %d cells with a spatial bin size of %dcm\n',nCells,params.SpatialBin);
        %%
%         allRho = nan(numel(cells_to_plot),1);
        l=1;
        plotWidth = 160*l;
        plotHeight = 500*l;
        h = figure('Position',[100 100 plotWidth plotHeight]); hold on; % changed from 160 to 320 for repeating, also from 500 to 166 for 50 trials
        for k = 1:numel(cells_to_plot)
            
            fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));
    
            % get spike times and index into post for cell k 
            spike_t = sp.st(sp.clu==cells_to_plot(k));
            [~,~,spike_idx] = histcounts(spike_t,post);
            spike_idx(spike_idx==0) = []; %remove spike indexes that don't exist in post (e.g. spikes that happen while the animal is stationary when post was speed filtered)
               
            singleCellallTrialsFR = nan(max(trial),spatialBins);
            for i = 1:max(trial)
                itrial_kfr = calcSmoothedFR_SpatialBin(spike_idx(i==trial(spike_idx)), posx(i==trial),posx, p, trackEnd);
                singleCellallTrialsFR(i,:) = itrial_kfr;
            end
            
            clf;
            imagesc(singleCellallTrialsFR)
            set(gca,'YDir','normal')
            colormap(hot(8))
            title(sprintf('c%d, d=%d',cells_to_plot(k),round(spike_depth(k))));
             % save fig
            if save_figs
                saveas(h,fullfile(image_save_dir,sprintf('%d.png',k)),'png');
                saveas(h,fullfile(image_save_dir,sprintf('%d.pdf',k)),'pdf');
            end

            % Calculate Stability of 2nd 3rd of baseline compared
            % last 3rd of baseline
%             fr1 = calcSmoothedFR_SpatialBin(spike_idx(trial(spike_idx)<40 & trial(spike_idx)>30), posx(trial<40 & trial>30),posx, p, trackEnd);
%             fr2 = calcSmoothedFR_SpatialBin(spike_idx(trial(spike_idx)<50 & trial(spike_idx)>40), posx(trial<50 & trial>40),posx, p, trackEnd);

%             fr1 = calcSmoothedFR_SpatialBin(spike_idx(trial(spike_idx)<preDrugTrials/2), posx(trial<preDrugTrials/2),posx, p, trackEnd);
%             fr2 = calcSmoothedFR_SpatialBin(spike_idx(trial(spike_idx)>preDrugTrials/2), posx(trial>preDrugTrials/2),posx, p, trackEnd);
%             rho = corr(fr1',fr2');
            
            testTrials = 1:50;
            rho = trialByTrialStability(testTrials, trial, spike_idx, posx, p, trackEnd);
            allRho(size(allRho,2)+1) = rho;
   
             %count total num of cells
            totalCell = totalCell + 1;
            
            % plot raster plot
            if rho > 0.01
                cla;


                if preDrugTrials ~= 0
                    if rho > 0.5
                        plot(posx(spike_idx),trial(spike_idx)-preDrugTrials,'r.');
                    else
                        plot(posx(spike_idx),trial(spike_idx)-preDrugTrials,'k.');
                    end
                    ylim([-preDrugTrials max(trial)+1-preDrugTrials]);
                else
                    plot(posx(spike_idx),trial(spike_idx),'k.');
                    ylim([0 max(trial)+1]);
                end
                
               
                xlim([params.TrackStart trackLength]);
                title(sprintf('c%d, d=%d, rho=%.3f\n %s %s',cells_to_plot(k),round(spike_depth(k)), rho,animalName, sessionDate ));
%                 xticks(''); yticks('');
                
                
                % save fig
                if save_figs
                    saveas(h,fullfile(image_save_dir,sprintf('%s_%s_%d.png',animalName,sessionDate,k)),'png');
                end
            end


        end
        

    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end
fprintf(strcat('\nTotal Number of Cells: ',num2str(totalCell),'\n'));
