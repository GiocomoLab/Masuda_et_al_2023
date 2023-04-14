%% LFP analysis
% First, run trans_prob_analysis to generate wt_ket_sessions.mat file and calculate_lfp_power.ipynb to generate lfp freq band csvs
load('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\wt_ket_sessions.mat')

[peri_ctrl, peri_ket, delta_ctrl, delta_ket] = deal([]);
output_dir = '//oak-smb-giocomo.stanford.edu/groups/giocomo/fkmasuda/fkm_analysis/EAJ_revisions/lfp/'
for s = 1:length(wt_ket_sessions)
    theta_power = load_table(output_dir+wt_ket_sessions{s,2}+'_theta_power.csv')
	
	% find 1s bins closest to injection times
    bins = 0:1:length(theta_power);
    [~, base_end] = min(abs(bins-wt_ket_sessions{s,3}));
    [~, ctrl_5_before] = min(abs(bins-(wt_ket_sessions{s,3}-5*60)));
    [~, ctrl_inj] = min(abs(bins-wt_ket_sessions{s,3}));
    [~, ctrl_5_after] = min(abs(bins-(wt_ket_sessions{s,3}+5*60)));
    [~, ctrl_10_after] = min(abs(bins-(wt_ket_sessions{s,3}+10*60)));
    [~, ket_5_before] = min(abs(bins-(wt_ket_sessions{s,4}-5*60)));
    [~, ket_inj] = min(abs(bins-wt_ket_sessions{s,4}));
    [~, ket_5_after] = min(abs(bins-(wt_ket_sessions{s,4}+5*60)));
    [~, ket_10_after] = min(abs(bins-(wt_ket_sessions{s,4}+10*60)));
	
	% get theta power in time bins
	peri_ctrl = [peri_ctrl; theta_power(ctrl_5_before:ctrl_10_after,:)'];
	peri_ket = [peri_ket; theta_power(ket_5_before:ket_10_after,:)'];
	delta_ctrl = [delta_ctrl; [nanmean(theta_power(ctrl_inj:ctrl_5_after,:)-...
        theta_power(ctrl_5_before:ctrl_inj,:))./nanmean(theta_power(1:base_end,:))]'];
	delta_ket = [delta_ket; [nanmean(theta_power(ket_inj:ket_5_after,:)-...
        theta_power(ket_5_before:ket_inj,:))./nanmean(theta_power(1:base_end,:))]'];
		
save('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\theta_power.mat')

%% Plot
% 3E/F
x = -5:.333:10;
plotBinnedDataOverTime(x, peri_ctrl, peri_ket 'Theta Power')

% 3G
figure; clf; clear g;
delta_ctrl(~isfinite(delta_ctrl)) = NaN;
delta_ket(~isfinite(delta_ket)) = NaN;
y = vertcat(delta_ctrl,delta_ket,);
color = vertcat(repmat({'Ctrl'},size(delta_ctrl_fromE1,1),1),...
    repmat({'Ket'},size(delta_ket_fromE1,1),1));
g=gramm('x',color,'y',y,'color',color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',3);
customColorMap = [
    0.8 0.2 0.8
    0 0.8 0.2];
g.set_color_options('map',customColorMap);
g.draw;

%% Statistics
compare_centers(delta_ctrl, delta_ket, 'paired', true)