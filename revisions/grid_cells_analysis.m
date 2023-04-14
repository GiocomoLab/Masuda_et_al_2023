addpath(genpath("C:\Users\Niflheim\Documents\GitHub\Giocomo\JohnKeiNPAnalysis"))
addpath(genpath('//oak-smb-giocomo.stanford.edu/groups/giocomo/export/data/Users/KMasuda/Neuropixels/plottingFxns'))
addpath(genpath("C:/Users/Niflheim/Documents/MATLAB/gramm"))

%% Repeating analyses with grid cells
%% Filter & save grid cells
% Run masterScript first, then...
% Filter for WT mec cells during ketamine sessions
seshIndx = ismember(allCells.metadata(:,8),'ketamine');
ketamineCells = filterAllCellsStruct(allCells,seshIndx);
seshIndx = ismember(ketamineCells.metadata(:,4),'WT');
wt_ket_cells = filterAllCellsStruct(ketamineCells,seshIndx);

% Filter for excitatory (<15Hz), gain change neurons
seshIndx = ~wt_ket_cells.interneuronFlag;
wt_ket_cells_noInterneurons = filterAllCellsStruct(wt_ket_cells,seshIndx);
seshIndx = logical(wt_ket_cells_noInterneurons.gainModulationValues(:,4));
wt_ket_cells_stableGainChange = filterAllCellsStruct(wt_ket_cells_noInterneurons,seshIndx);

%% Additional spatial stability & spatial information filters
stability = zeros(length(wt_ket_cells_stableGainChange.metadata),1);
for c = 1:length(wt_ket_cells_stableGainChange.metadata)
    baseline_stability = wt_ket_cells_stableGainChange.stabilityScoreCurve(c,1:50);
    baseline_stability(~isfinite(baseline_stability)) = 0;
    stability(c) = mean(baseline_stability);
end
seshIndx = stability>0.2;
grid_cells = filterAllCellsStruct(wt_ket_cells_stableGainChange, seshIndx);
baseline_SI = mean(grid_cells.bitsPerSecCurve(:,1:50),2);
seshIndx = baseline_SI>3;
grid_cells = filterAllCellsStruct(grid_cells,seshIndx);
save('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\grid_cells.mat', 'grid_cells', '-v7.3')

%% Plot
% 3I-J-K
plot_peakinessCurves(grid_cells)
% 3M
plot_timeWarpedStabilityScoreCurves(grid_cells)
% 5C
plot_AttractorDynamicsbySession_figure(grid_cells)

%% Grid Scale: shuffle method
smoothSigma = 2;
smoothWindow = floor(smoothSigma*5/2)*2+1;
gauss_filter = fspecial('gaussian',[smoothWindow 1], smoothSigma);

trial_blocks = 1:50:290;
trial_blocks = [trial_blocks, 290];
n_cells = length(grid_cells.metadata);
[grid_scale, n_fields, field_width] = deal(zeros(n_cells, length(trial_blocks)-1));
for t = 1:length(trial_blocks)-1
    for c = 1:n_cells %iterate over all cells
        spatial_map = mean(squeeze(grid_cells.spatialFRsmooth(c,trial_blocks(t):trial_blocks(t+1),:)));
        shuffle_prominence = [];
        shuffle_peak = [];
        % find peaks of shuffled distribution
        for s = 1:1000
            % shuffle, then smooth
            unsmoothed_map = mean(squeeze(grid_cells.spatialFR2(c,trial_blocks(t):trial_blocks(t+1),:)));
            shuffle = unsmoothed_map(randperm(length(unsmoothed_map)));
            % smooth & find peaks over repeated map because the track is
            % circular (infinite loop)
            shuffle_smooth = conv(repmat(shuffle,1,3),gauss_filter,'same');
            [peak,bins,~,prominence] = findpeaks(shuffle_smooth, 'MinPeakWidth', 3);
            % subset on peaks within non-repeated (middle copy) data
            idx = bins>length(shuffle) & bins<length(shuffle)*2+1;
            shuffle_prominence = [shuffle_prominence, prominence(idx)];
            shuffle_peak = [shuffle_peak, peak(idx)];
        end
        if ~isempty(shuffle_prominence)
            min_prominence = prctile(shuffle_prominence, 80); %greater prominence than X% of shuffled peaks
            min_peak = prctile(shuffle_peak, 80); %greater prominence than X% of shuffled peaks
        else
            min_prominence = 0;
            min_peak = 0;
        end
        [~, peak_bins, widths, ~] = findpeaks(repmat(spatial_map,1,3), 'MinPeakHeight', min_peak, 'MinPeakWidth', 3);
        if max(widths)>100
            figure
            findpeaks(repmat(spatial_map,1,3), 'MinPeakProminence', min_prominence, 'MinPeakWidth', 3)
        end
%         [~, peak_bins, widths, ~] = findpeaks(repmat(spatial_map,1,3), 'MinPeakProminence', min_prominence, 'MinPeakWidth', 3);
%         if max(widths)>50
%             disp(max(widths))
%         end
        idx = peak_bins>length(spatial_map) & peak_bins<length(spatial_map)*2+1;
        n_fields(c,t) = sum(idx);
        if isempty(idx)
            field_width(c,t) = NaN;
        else
            field_width(c,t) = max(widths(idx));
        end
        if n_fields(c,t)>1
            grid_scale(c,t) = min(peak_bins(2:end)-peak_bins(1:end-1))*2;
        else
            grid_scale(c,t) = 400;
        end
    end
end

save('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\grid_scale.mat', 'grid_scale')
save('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\n_fields.mat', 'n_fields')
save('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\field_width.mat', 'field_width')

%% Grid Scale: FR percentile method
% % too strict with height requirement of 3Hz, too permissive without
% trial_blocks = 1:50:290;
% trial_blocks = [trial_blocks, 290];
% n_cells = length(grid_cells.metadata);
% grid_scale = zeros(n_cells, length(trial_blocks)-1);
% n_fields = zeros(n_cells, length(trial_blocks)-1);
% for t = 1:length(trial_blocks)-1
%     for c = 1%:n_cells %iterate over all cells
%         fr = squeeze(grid_cells.spatialFRsmooth(c,trial_blocks(t):trial_blocks(t+1),:));
%         spatial_map = mean(fr);
%         min_prominence = 0.05*prctile(fr(:), 99.5);
%         [~, peak_bins] = findpeaks(spatial_map, 'MinPeakWidth', 5, 'MinPeakProminence', min_prominence);
%         figure
%         findpeaks(spatial_map, 'MinPeakWidth', 5, 'MinPeakProminence', min_prominence)
%         n_fields(c,t) = length(peak_bins);
%         if n_fields(c,t)>1
%             grid_scale(c,t) = min(peak_bins(2:end)-peak_bins(1:end-1))*2; %in cm
%         else
%             grid_scale(c,t) = 400;
%         end
%     end
% end

%% Plot
figure; clf; clear g;
y = grid_scale(:);
color = vertcat(repmat({'Baseline'},n_cells,1), repmat({'Control'},n_cells,1),...
    repmat({'Early Ketamine'},n_cells,1), repmat({'Early Mid Ketamine'},n_cells,1),...
    repmat({'Late Mid Ketamine'},n_cells,1), repmat({'Late Ketamine'},n_cells,1));
g=gramm('x',color,'y',y,'color',color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',0.8);
g.draw;
% 
figure; clf; clear g;
y = n_fields(:);
g=gramm('x',color,'y',y,'color',color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',0.8);
g.draw;

figure; clf; clear g;
y = field_width(:);
g=gramm('x',color,'y',y,'color',color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',0.8);
g.draw;

figure; clf; clear g;
control_scale = grid_scale(:,2)-grid_scale(:,1);
control_scale = control_scale(n_fields(:,1)>0 & n_fields(:,2)>0);
early_scale = grid_scale(:,3)-grid_scale(:,1);
early_scale = early_scale(n_fields(:,1)>0 & n_fields(:,3)>0);
% earlymid_scale = grid_scale(:,4)-grid_scale(:,1);
% earlymid_scale = earlymid_scale(n_fields(:,1)>0 & n_fields(:,4)>0);
% midlate_scale = grid_scale(:,5)-grid_scale(:,1);
% midlate_scale = midlate_scale(n_fields(:,1)>0 & n_fields(:,5)>0);
verylate_scale = grid_scale(:,6)-grid_scale(:,1);
verylate_scale = verylate_scale(n_fields(:,1)>0 & n_fields(:,6)>0);
% y = vertcat(control_scale, early_scale, earlymid_scale, midlate_scale, verylate_scale);
y = vertcat(control_scale, early_scale, verylate_scale);
% color = vertcat(repmat({'Control'},length(control_scale),1), ...
%     repmat({'Early Ketamine'},length(early_scale),1), repmat({'Early Mid Ketamine'},length(earlymid_scale),1),...
%     repmat({'Late Mid Ketamine'},length(midlate_scale),1), repmat({'Very Late Ketamine'},length(verylate_scale),1));
color = vertcat(repmat({'Control'},length(control_scale),1), ...
    repmat({'Early Ketamine'},length(early_scale),1),...
    repmat({'Very Late Ketamine'},length(verylate_scale),1));
g=gramm('x',color,'y',y,'color',color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',0.8);
g.draw;

figure; clf; clear g;
control_fields = n_fields(:,2)-n_fields(:,1);
control_fields = control_fields(n_fields(:,1)>0 & n_fields(:,2)>0);
early_fields = n_fields(:,3)-n_fields(:,1);
early_fields = early_fields(n_fields(:,1)>0 & n_fields(:,3)>0);
% earlymid_fields = n_fields(:,4)-n_fields(:,1);
% midlate_fields = n_fields(:,5)-n_fields(:,1);
verylate_fields = n_fields(:,6)-n_fields(:,1);
verylate_fields = verylate_fields(n_fields(:,1)>0 & n_fields(:,6)>0);
y = vertcat(control_fields, early_fields, verylate_fields);
g=gramm('x',color,'y',y,'color',color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',0.8);
g.draw;

figure; clf; clear g;
control_width = field_width(:,2)-field_width(:,1);
control_width = control_width(n_fields(:,1)>0 & n_fields(:,2)>0);
early_width = field_width(:,3)-field_width(:,1);
early_width = early_width(n_fields(:,1)>0 & n_fields(:,3)>0);
% earlymid_width = field_width(:,4)-field_width(:,1);
% earlymid_width = earlymid_width(n_fields(:,1)>0 & n_fields(:,4)>0);
% midlate_width = field_width(:,5)-field_width(:,1);
% midlate_width = midlate_width(n_fields(:,1)>0 & n_fields(:,5)>0);
verylate_width = field_width(:,6)-field_width(:,1);
verylate_width = verylate_width(n_fields(:,1)>0 & n_fields(:,6)>0);
y = vertcat(control_width, early_width, verylate_width);
% y = vertcat(control_width, early_width, earlymid_width, midlate_width, verylate_width);
g=gramm('x',color,'y',y,'color',color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',0.8);
g.draw;

%% Plot example ratemaps
figure
hold on
for t = 1:length(trial_blocks)-1
    for c = 55
        spatial_map = mean(squeeze(grid_cells.spatialFRsmooth(c,trial_blocks(t):trial_blocks(t+1),:)));
        shuffle_prominence = [];
        shuffle_peak = [];
        % find peaks of shuffled distribution
        for s = 1:1000
            % shuffle, then smooth
            unsmoothed_map = mean(squeeze(grid_cells.spatialFR2(c,trial_blocks(t):trial_blocks(t+1),:)));
            shuffle = unsmoothed_map(randperm(length(unsmoothed_map)));
            % smooth & find peaks over repeated map because the track is
            % circular (infinite loop)
            shuffle_smooth = conv(repmat(shuffle,1,3),gauss_filter,'same');
            [peak,bins,~,prominence] = findpeaks(shuffle_smooth, 'MinPeakWidth', 3);
            % subset on peaks within non-repeated (middle copy) data
            idx = bins>length(shuffle) & bins<length(shuffle)*2+1;
            shuffle_prominence = [shuffle_prominence, prominence(idx)];
            shuffle_peak = [shuffle_peak, peak(idx)];
        end
        if ~isempty(shuffle_prominence)
            min_prominence = prctile(shuffle_prominence, 80); %greater prominence than X% of shuffled peaks
            min_peak = prctile(shuffle_peak, 80); %greater prominence than X% of shuffled peaks
        else
            min_prominence = 0;
            min_peak = 0;
        end
        subplot(6,1,t);
        findpeaks(repmat(spatial_map,1,3), 'MinPeakProminence', min_prominence, 'MinPeakWidth', 3)
%         findpeaks(repmat(spatial_map,1,3), 'MinPeakHeight', min_peak, 'MinPeakWidth', 3)
        [~, peak_bins, widths, ~] = findpeaks(repmat(spatial_map,1,3), 'MinPeakProminence', min_prominence, 'MinPeakWidth', 3);
        xlim([200 400])
        idx = peak_bins>length(spatial_map) & peak_bins<length(spatial_map)*2+1;
        disp(mean(widths(idx)))
    end
end
hold off

%% Plot example rasters
image_dir = '//oak-smb-giocomo.stanford.edu/groups/giocomo/fkmasuda/fkm_analysis/EAJ_revisions/rasters';
for c = 1:n_cells
    h = figure('Position',[100 100 160 500]); hold on;
    posx = grid_cells.posX(c);
    trial = grid_cells.trial(c);
    spike_idx = grid_cells.spike_idx{c};
    plot(posx.posx(spike_idx),trial.trial(spike_idx),'k.', 'MarkerSize', 4)
    xlim([0 400])
    ylim([0 290])
    set(gca, 'YDir', 'reverse')
    title(c);
    saveas(h,fullfile(image_dir,sprintf('%d.png',c)),'png');
end
combineRASTERS("All_grid_cells", 781, 250, 6, image_dir)

%% Statistics
% Wilcoxon matched pairs signed rank test
% compare_centers(control_width, early_width)
compare_centers(control_width(n_fields(n_fields(:,1)>0 & n_fields(:,2)>0,3)>0),...
    early_width(n_fields(n_fields(:,1)>0 & n_fields(:,3)>0,2)>0), 'paired', true, 'multcompare', 3)
% compare_centers(early_width, verylate_width)
compare_centers(early_width(n_fields(n_fields(:,1)>0 & n_fields(:,3)>0,6)>0),...
    verylate_width(n_fields(n_fields(:,1)>0 & n_fields(:,6)>0,3)>0), 'paired', true, 'multcompare', 3)
% compare_centers(control_width, verylate_width)
compare_centers(control_width(n_fields(n_fields(:,1)>0 & n_fields(:,2)>0,6)>0),...
    verylate_width(n_fields(n_fields(:,1)>0 & n_fields(:,6)>0,2)>0), 'paired', true, 'multcompare', 3)