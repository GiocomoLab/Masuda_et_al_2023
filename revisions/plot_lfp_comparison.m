%% LFP analysis
% First, calculate_lfp_power.ipynb to generate lfp freq band csvs
load('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\wt_ket_sessions.mat')

[peri_ctrl, peri_ket, delta_ctrl, delta_ket] = deal([]);
[base_trace, ctrl_trace, ket_trace] = deal(NaN(length(wt_ket_sessions), 15000));
output_dir = '//oak-smb-giocomo.stanford.edu/groups/giocomo/fkmasuda/fkm_analysis/EAJ_revisions/lfp/';
for s = 1:length(wt_ket_sessions)
    theta_power = table2array(readtable([output_dir,wt_ket_sessions{s,2},'_theta_power.csv']));
	
	% find 1s bins closest to injection times
    bins = 1:length(theta_power);
    [~, base_end] = min(abs(bins-wt_ket_sessions{s,3}));
    [~, ctrl_end] = min(abs(bins-wt_ket_sessions{s,4}));
    [~, ctrl_5_before] = min(abs(bins-(wt_ket_sessions{s,3}-5*60)));
    [~, ctrl_inj] = min(abs(bins-wt_ket_sessions{s,3}));
    [~, ctrl_5_after] = min(abs(bins-(wt_ket_sessions{s,3}+5*60)));
    [~, ctrl_10_after] = min(abs(bins-(wt_ket_sessions{s,3}+10*60)));
    [~, ket_5_before] = min(abs(bins-(wt_ket_sessions{s,4}-5*60)));
    [~, ket_inj] = min(abs(bins-wt_ket_sessions{s,4}));
    [~, ket_5_after] = min(abs(bins-(wt_ket_sessions{s,4}+5*60)));
    [~, ket_10_after] = min(abs(bins-(wt_ket_sessions{s,4}+10*60)));

    % normalize
    theta_power_norm = ((theta_power-nanmean(theta_power(1:base_end,:)))/nanstd(theta_power(1:base_end,:)));
	
	% get theta power in time bins
	peri_ctrl = [peri_ctrl; theta_power_norm(ctrl_5_before:ctrl_10_after,:)'];
	peri_ket = [peri_ket; theta_power_norm(ket_5_before:ket_10_after,:)'];
	delta_ctrl = [delta_ctrl; nanmean(theta_power_norm(ctrl_inj:ctrl_5_after,:)-...
        nanmean(theta_power_norm(ctrl_5_before:ctrl_inj,:)))];
	delta_ket = [delta_ket; nanmean(theta_power_norm(ket_inj:ket_5_after,:)-...
        nanmean(theta_power_norm(ket_5_before:ket_inj,:)))];

    base_trace(s,1:length(theta_power_norm(1:base_end,:))) = theta_power_norm(1:base_end,:)';
    ctrl_trace(s,1:length(theta_power_norm(base_end+1:ctrl_end,:))) = theta_power_norm(base_end+1:ctrl_end,:)';
    ket_trace(s,1:length(theta_power_norm(ctrl_end+1:end,:))) = theta_power_norm(ctrl_end+1:end,:)';

    % get velocity
    data_file = ['//oak-smb-giocomo.stanford.edu/groups/giocomo/fkmasuda/fkm_analysis/combinedSesh/fr_data_matrices_noSmoothing/',...
        wt_ket_sessions{s,2},'_baseline1+controlinjx1+ketamine1_fr+corr.mat'];
    load(data_file,'speed')
%     disp((length(speed)/50-length(theta_power))/length(theta_power))
end
%save('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\theta_power.mat')

%% Plot and statistics
% 3E/F
% x = -5:1/60:10;
% plotBinnedDataOverTime(x, peri_ctrl, peri_ket)
% x = 1/60:1/60:30;
% plotBinnedDataOverTime(x, ctrl_trace(:,1:length(x)), ket_trace(:,1:length(x)), 'base', base_trace(:,1:length(x)))

rescale = 60;
x = 1:30;
for a = 1:size(ket,1)
    base_smoothed(a,:) = nanmean(reshape(base_trace(a,1:rescale*fix(length(base_trace)/rescale)),rescale,[]));
    ctrl_smoothed(a,:) = nanmean(reshape(ctrl_trace(a,1:rescale*fix(length(ctrl_trace)/rescale)),rescale,[]));
    ket_smoothed(a,:) = nanmean(reshape(ket_trace(a,1:rescale*fix(length(ket_trace)/rescale)),rescale,[]));
end
plotBinnedDataOverTime(x, ctrl_smoothed(:,1:length(x)), ket_smoothed(:,1:length(x)), 'base', base_smoothed(:,1:length(x)))

x = 1/60:1/60:30;
mean_base = nanmean(base_trace(:,1:length(x)),2);
mean_ctrl = nanmean(ctrl_trace(:,1:length(x)),2);
mean_ket = nanmean(ket_trace(:,1:length(x)),2);

% 3G
% figure; clf; clear g;
% y = vertcat(mean_base, mean_ctrl, mean_ket);
% color = vertcat(repmat({'Base'},size(mean_base,1),1),...
%     repmat({'Ctrl'},size(mean_ctrl,1),1),...
%     repmat({'Ket'},size(mean_ket,1),1));
% g=gramm('x',color,'y',y,'color',color);
% g.stat_violin('normalization','width','dodge',0,'fill','edge');
% g.stat_boxplot('width',1);
% customColorMap = [ 0.5 0.5 0.5
%     0.8 0.2 0.8
%     0 0.8 0.2];
% g.set_color_options('map',customColorMap);
% g.draw;

% compare_centers(delta_ctrl, delta_ket, 'paired', true)
compare_centers(mean_base, mean_ctrl, 'paired', true)
compare_centers(mean_ctrl, mean_ket, 'paired', true)
compare_centers(mean_base, mean_ket, 'paired', true)