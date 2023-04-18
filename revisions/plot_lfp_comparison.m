%% LFP analysis
% First, calculate_lfp_power.ipynb to generate lfp freq band csvs
load('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\wt_ket_sessions.mat')

[peri_ctrl, peri_ket, delta_ctrl, delta_ket] = deal([]);
[base_trace, ctrl_trace, ket_trace] = deal(NaN(length(wt_ket_sessions), 15000));
output_dir = '//oak-smb-giocomo.stanford.edu/groups/giocomo/fkmasuda/fkm_analysis/EAJ_revisions/lfp/';
for s = 1:length(wt_ket_sessions)
    theta_power = table2array(readtable([output_dir,wt_ket_sessions{s,2},'_FG_power.csv']));

    % find 1s bins closest to injection times
    bins = 1:length(theta_power);
    [~, base_start] = min(abs(bins-wt_ket_sessions{s,5}));
    [~, base_end] = min(abs(bins-wt_ket_sessions{s,6}));
    [~, ctrl_start] = min(abs(bins-wt_ket_sessions{s,7}));
    [~, ctrl_end] = min(abs(bins-wt_ket_sessions{s,8}));
    [~, ket_start] = min(abs(bins-wt_ket_sessions{s,9}));
    [~, ket_end] = min(abs(bins-wt_ket_sessions{s,10}));

    [~, ctrl_5_before] = min(abs(bins-(wt_ket_sessions{s,7}-5*60)));
    [~, ctrl_inj] = min(abs(bins-wt_ket_sessions{s,7}));
    [~, ctrl_5_after] = min(abs(bins-(wt_ket_sessions{s,7}+5*60)));
    [~, ctrl_10_after] = min(abs(bins-(wt_ket_sessions{s,7}+10*60)));
    [~, ket_5_before] = min(abs(bins-(wt_ket_sessions{s,9}-5*60)));
    [~, ket_inj] = min(abs(bins-wt_ket_sessions{s,9}));
    [~, ket_5_after] = min(abs(bins-(wt_ket_sessions{s,9}+5*60)));
    [~, ket_10_after] = min(abs(bins-(wt_ket_sessions{s,9}+10*60)));

    % normalize
    theta_power_norm = ((theta_power-nanmean(theta_power(base_start:base_end,:)))/...
        nanstd(theta_power(base_start:base_end,:)));

    % get theta power in time bins
    peri_ctrl = [peri_ctrl; theta_power_norm(ctrl_5_before:ctrl_10_after,:)'];
    peri_ket = [peri_ket; theta_power_norm(ket_5_before:ket_10_after,:)'];
    delta_ctrl = [delta_ctrl; nanmean(theta_power_norm(ctrl_inj:ctrl_5_after,:)-...
        nanmean(theta_power_norm(ctrl_5_before:ctrl_inj,:)))];
    delta_ket = [delta_ket; nanmean(theta_power_norm(ket_inj:ket_5_after,:)-...
        nanmean(theta_power_norm(ket_5_before:ket_inj,:)))];

    base_trace(s,1:length(theta_power_norm(base_start:base_end,:))) = theta_power_norm(base_start:base_end,:)';
    ctrl_trace(s,1:length(theta_power_norm(ctrl_start:ctrl_end,:))) = theta_power_norm(ctrl_start:ctrl_end,:)';
    ket_trace(s,1:length(theta_power_norm(ket_start:ket_end,:))) = theta_power_norm(ket_start:ket_end,:)';
end
%save('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\theta_power.mat')

%% Plot and statistics
% 3E/F
rescale = 60;
% x = -5:9;
% for a = 1:length(wt_ket_sessions)
%     peri_ctrl_smoothed(a,:) = nanmean(reshape(peri_ctrl(a,1:rescale*fix(length(peri_ctrl)/rescale)),rescale,[]));
%     peri_ket_smoothed(a,:) = nanmean(reshape(peri_ket(a,1:rescale*fix(length(peri_ket)/rescale)),rescale,[]));
% end
% plotBinnedDataOverTime(x, peri_ctrl_smoothed, peri_ket_smoothed)

x = 1:30;
for a = 1:length(wt_ket_sessions)
    base_smoothed(a,:) = nanmean(reshape(base_trace(a,1:rescale*fix(length(base_trace)/rescale)),rescale,[]));
    ctrl_smoothed(a,:) = nanmean(reshape(ctrl_trace(a,1:rescale*fix(length(ctrl_trace)/rescale)),rescale,[]));
    ket_smoothed(a,:) = nanmean(reshape(ket_trace(a,1:rescale*fix(length(ket_trace)/rescale)),rescale,[]));
end
plotBinnedDataOverTime(x, ctrl_smoothed(:,1:length(x)), ket_smoothed(:,1:length(x)), 'base', base_smoothed(:,1:length(x)))

%% Match based on velocity
x = 1/60:1/60:30;
[base_trace_speed, ctrl_trace_speed, ket_trace_speed] = deal(NaN(length(wt_ket_sessions), length(x)));
for s = 1:length(wt_ket_sessions)
    % get velocity
    data_file = ['//oak-smb-giocomo.stanford.edu/groups/giocomo/fkmasuda/fkm_analysis/combinedSesh/fr_data_matrices_noSmoothing/',...
        wt_ket_sessions{s,2},'_baseline1+controlinjx1+ketamine1_fr+corr.mat'];
    load(data_file,'speed')

    %change from 20ms speed bins to 1s LFP bins
    rescale = 50;
    speed = nanmean(reshape(speed(1:rescale*fix(length(speed)/rescale)),rescale,[]));

    % get indices of speed vector (different from indices of lfp)
    ctrl_idx = int16(wt_ket_sessions{s,3});
    ket_idx = int16(wt_ket_sessions{s,4});
    base_end = min(ctrl_idx, 1800);
    ctrl_end = min(ket_idx, ctrl_idx+1800-1);
    ket_end = min(length(speed), ket_idx+1800-1);
    base_speed = speed(1:base_end);
    ctrl_speed = speed(ctrl_idx:ctrl_end);
    ket_speed = speed(ket_idx:ket_end);

    % sample ctrl and ket speed vectors so they match base speed vectors
    [unique_ctrl_speed, return_idx, ~] = unique(ctrl_speed, 'stable');
    ctrl_match_base_idx = interp1(unique_ctrl_speed, 1:length(unique_ctrl_speed), base_speed, 'nearest', 'extrap');
    ctrl_match_base_idx = return_idx(ctrl_match_base_idx);
    [unique_ket_speed, return_idx, ~] = unique(ket_speed, 'stable');
    ket_match_base_idx = interp1(unique_ket_speed, 1:length(unique_ket_speed), base_speed, 'nearest', 'extrap');
    ket_match_base_idx = return_idx(ket_match_base_idx);

    % get ctrl and ket lfp samples so speed is matched to base
    base_trace_speed(s,1:length(base_speed)) = base_trace(s,1:length(base_speed));
    ctrl_trace_speed(s,1:length(ctrl_match_base_idx)) = ctrl_trace(s,ctrl_match_base_idx);
    ket_trace_speed(s,1:length(ket_match_base_idx)) = ket_trace(s,ket_match_base_idx);
end

mean_base_speed = nanmean(base_trace_speed,2);
mean_ctrl_speed = nanmean(ctrl_trace_speed,2);
mean_ket_speed = nanmean(ket_trace_speed,2);

figure; clf; clear g;
y = vertcat(mean_base_speed, mean_ctrl_speed, mean_ket_speed);
color = vertcat(repmat({'Base'},size(mean_base_speed,1),1),...
    repmat({'Ctrl'},size(mean_ctrl_speed,1),1),...
    repmat({'Ket'},size(mean_ket_speed,1),1));
g=gramm('x',color,'y',y,'color',color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',1);
customColorMap = [ 0.5 0.5 0.5
    0.8 0.2 0.8
    0 0.8 0.2];
g.set_color_options('map',customColorMap);
g.draw;

stats = compare_centers(mean_base_speed, mean_ctrl_speed, 'paired', true);
fprintf("Base vs Ctrl: %s\n",stats.resultsstr)
stats = compare_centers(mean_ctrl_speed, mean_ket_speed, 'paired', true);
fprintf("Ctrl vs Ket: %s\n",stats.resultsstr)
stats = compare_centers(mean_base_speed, mean_ket_speed, 'paired', true);
fprintf("Base vs Ket: %s\n",stats.resultsstr)