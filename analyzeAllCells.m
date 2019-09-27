        
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));

addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/functions'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/spikes'));
% some params
params = readtable('UniversalParams.xlsx');
p = params;
save_figs = false;
%%

% sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/*.mat');
sessions = dir('/Users/KeiMasuda/Desktop/fkm_analysis/fr_corr_matrices_noSpeedFilter/*.mat');
% Get subset of desired sessions
sessions = filterSessions(sessions, 'WT');

%%
totalCell = 0;
totalSpatialCell = 0;
allRho = [];
allSpatialIndx = {};
for n = 1:numel(sessions)
    try
        
        matPath = fullfile(sessions(n).folder, sessions(n).name);
        session_name = sessions(n).name(1:end-4);
        animalName = extractBefore(session_name,'_');
        sessionDate = extractBefore(extractAfter(session_name,'_'),'_');
        trackLength = 400;
        load(fullfile(matPath), 'all_fr', 'avg_all_fr', 'all_corrmatrix', 'avg_all_corrmatrix', ...
             'all_waveforms', 'cells_to_plot','spike_depth','all_drugEffectScores',...
            'trial','all_cellCorrScore','trials_corrTemplate', 'avg_all_cellCorrScore', 'avg_cell_fr');
        
%         imgDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/img';
%         imgDir = '/Users/KeiMasuda/Desktop/fkm_analysis/img/ratemaps';
        imgDir = '/Users/KeiMasuda/Desktop/fkm_analysis/img/rhoGreater01';
        nCells = numel(cells_to_plot);
        totalCell = totalCell + nCells;
        fprintf('Calculating firing rate for %d cells with a spatial bin size of %dcm\n',nCells,params.SpatialBin);
        %%
        close all;
        spatialIndx = [];
        l=1;
        plotWidth = 160*l;
        plotHeight = 500*l;
        h = figure('Position',[100 100 plotWidth plotHeight]); hold on; % changed from 160 to 320 for repeating, also from 500 to 166 for 50 trials
        for k = 1:numel(cells_to_plot)
            
            
            frMap = squeeze(all_fr(k,:,:));
            % Calculate Stability of 
            testTrials = 1:50;
            numAvgTrials = 10; % compare 1:5 with 6:10 for ever 10 testTrials
            rho = trialByTrialStability(frMap, testTrials, numAvgTrials);
            allRho(size(allRho,2)+1) = rho;
   
            clf; cla;
            imagesc(frMap)
            set(gca,'YDir','normal')
%             
            title(sprintf('%s-%s,\nc%d, d=%d \n rho=%.3f',animalName,sessionDate,cells_to_plot(k),round(spike_depth(k)),rho));
            
            % plot raster plot
            rhoThresh = 0.1;
            % rho > rhoThresh
            if rho > rhoThresh
                colormap(hot(10))
                totalSpatialCell = totalSpatialCell + 1;
                spatialIndx(size(spatialIndx,2)+1) = cells_to_plot(k);
                % save fig
                if save_figs
                    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));
                    saveas(h,fullfile(imgDir,sprintf('%s%s%s%s%d%s.png',animalName,'_',sessionDate,'_',k,'_spatial')),'png');
                end
            else
                colormap('default')
                if save_figs
                    fprintf('cell %d (%d/%d)\n',cells_to_plot(k),k,numel(cells_to_plot));
                    saveas(h,fullfile(imgDir,sprintf('%s%s%s%s%d.png',animalName,'_',sessionDate,'_',k)),'png');
                end
            end
            
        end
        seshIndx = sprintf('%s_%s',animalName,sessionDate);
        allSpatialIndx.(seshIndx) = spatialIndx;

    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end
fprintf(strcat('\nTotal Number of Cells: ',num2str(totalCell),'\nNumber of Spatial Cells:',num2str(totalSpatialCell),'\n'));
