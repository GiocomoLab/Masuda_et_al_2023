%% Synaptic transmission probability
%% Generate using CellExplorer (takes ~20mins per session)
% requires EAJ fork of CellExplorer: https://github.com/emilyasterjones/CellExplorer
run_CellExplorer_Masuda2023

%% Read results
load('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\wt_ket_sessions.mat')
[peri_ctrl_fromE, peri_ctrl_fromI, peri_ket_fromE, peri_ket_fromI,...
    delta_ctrl_fromE, delta_ctrl_fromI, delta_ket_fromE, delta_ket_fromI] = deal([]);
% [peri_ctrl_EE, peri_ctrl_EI, peri_ctrl_EU, peri_ctrl_IE, peri_ctrl_II, peri_ctrl_IU,...
%     peri_ctrl_fromE, peri_ctrl_toE, peri_ctrl_fromI, peri_ctrl_toI] = deal([]);
% [peri_ket_EE, peri_ket_EI, peri_ket_EU, peri_ket_IE, peri_ket_II, peri_ket_IU,...
%     peri_ket_fromE, peri_ket_toE, peri_ket_fromI, peri_ket_toI] = deal([]);
% [delta_ket_EE, delta_ctrl_EE, delta_ket_EI, delta_ctrl_EI, delta_ket_EU, delta_ctrl_EU, ...
%     delta_ket_IE, delta_ctrl_IE, delta_ket_II, delta_ctrl_II, delta_ket_IU, delta_ctrl_IU, ...
%     delta_ket_fromE, delta_ctrl_fromE, delta_ket_toE, delta_ctrl_toE, ...
%     delta_ket_fromI, delta_ctrl_fromI, delta_ket_toI, delta_ctrl_toI] = deal([]);
% [count_ee, count_ei, count_e, count_ie, count_ii, count_i] = deal([]);
[base_trace_fromE, ctrl_trace_fromE, ket_trace_fromE,...
    base_trace_fromI, ctrl_trace_fromI, ket_trace_fromI] = deal(NaN(10000, 180));
basepath = "\\oak-smb-giocomo.stanford.edu\groups\giocomo\export\data\Projects\JohnKei_NPH3";
count_e = 0;
count_i = 0;
for s = 1:length(wt_ket_sessions)
    metrics_file = sprintf('%s\\%s\\%s_keicontrasttrack_ketamine1_g0\\%s_keicontrasttrack_ketamine1_g0_imec0\\%s_keicontrasttrack_ketamine1_g0_imec0.cell_metrics.cellinfo',...
        basepath,wt_ket_sessions{s,1},wt_ket_sessions{s,2},wt_ket_sessions{s,2},wt_ket_sessions{s,2});
    load(metrics_file)

    % get indices of connections based on type (exc or inh or unk)
    fromE_connect = cell_metrics.putativeConnections.excitatory;
    fromI_connect = cell_metrics.putativeConnections.inhibitory;
%     [~, ~, exc_to_exc_idx] = intersect(fromE_connect(:,1), fromE_connect(:,2));
%     [~, ~, exc_to_inh_idx] = intersect(fromI_connect(:,1), fromE_connect(:,2));
%     exc_to_unk_idx = setxor(1:length(fromE_connect),[exc_to_exc_idx; exc_to_inh_idx]);
%     [~, ~, inh_to_exc_idx] = intersect(fromE_connect(:,1), fromI_connect(:,2));
%     [~, ~, inh_to_inh_idx] = intersect(fromI_connect(:,1), fromI_connect(:,2));
%     inh_to_unk_idx = setxor(1:length(fromI_connect),[inh_to_exc_idx; inh_to_inh_idx]);
%     count_e = [count_e; length(fromE_connect)];
%     count_i = [count_i; length(fromE_connect)];
%     count_ee = [count_ee; length(exc_to_exc_idx)];
%     count_ei = [count_ei; length(exc_to_inh_idx)];
%     count_ie = [count_ie; length(inh_to_exc_idx)];
%     count_ii = [count_ii; length(inh_to_inh_idx)];
    
    % get transmission probabilities by type
    fromE_trans_prob = cell_metrics.putativeConnections.excitatoryTransProbPriorCorr_binned;
%     EE_trans_prob = fromE_trans_prob(:,exc_to_exc_idx);
%     EI_trans_prob = fromE_trans_prob(:,exc_to_inh_idx);
%     EU_trans_prob = fromE_trans_prob(:,exc_to_unk_idx);
    fromI_trans_prob = cell_metrics.putativeConnections.inhibitoryTransProbPriorCorr_binned;
%     IE_trans_prob = fromI_trans_prob(:,inh_to_exc_idx);
%     II_trans_prob = fromI_trans_prob(:,inh_to_inh_idx);
%     IU_trans_prob = fromI_trans_prob(:,inh_to_unk_idx);

    %remove Inf, manually, because MATLAB is to stupid to do this
    fromE_trans_prob(~isfinite(fromE_trans_prob)) = NaN;
%     EE_trans_prob(~isfinite(EE_trans_prob)) = NaN;
%     EI_trans_prob(~isfinite(EI_trans_prob)) = NaN;
%     EU_trans_prob(~isfinite(EU_trans_prob)) = NaN;
    fromI_trans_prob(~isfinite(fromI_trans_prob)) = NaN;
%     IE_trans_prob(~isfinite(IE_trans_prob)) = NaN;
%     II_trans_prob(~isfinite(II_trans_prob)) = NaN;
%     IU_trans_prob(~isfinite(IU_trans_prob)) = NaN;

    % find 10s bins closest to injection times
    bins = cell_metrics.putativeConnections.bins;
    [~, base_start] = min(abs(bins-wt_ket_sessions{s,5}));
    [~, base_end] = min(abs(bins-wt_ket_sessions{s,6}));
    base_end = min(base_end, base_start+180);
    [~, ctrl_start] = min(abs(bins-wt_ket_sessions{s,7}));
    [~, ctrl_end] = min(abs(bins-wt_ket_sessions{s,8}));
    ctrl_end = min(ctrl_end, ctrl_start+180);
    [~, ket_start] = min(abs(bins-wt_ket_sessions{s,9}));
    [~, ket_end] = min(abs(bins-wt_ket_sessions{s,10}));
    ket_end = min(ket_end, ket_start+180);

    [~, ctrl_5_before] = min(abs(bins-(wt_ket_sessions{s,7}-5*60)));
    [~, ctrl_inj] = min(abs(bins-wt_ket_sessions{s,7}));
    [~, ctrl_5_after] = min(abs(bins-(wt_ket_sessions{s,7}+5*60)));
    [~, ctrl_10_after] = min(abs(bins-(wt_ket_sessions{s,7}+10*60)));
    [~, ket_5_before] = min(abs(bins-(wt_ket_sessions{s,9}-5*60)));
    [~, ket_inj] = min(abs(bins-wt_ket_sessions{s,9}));
    [~, ket_5_after] = min(abs(bins-(wt_ket_sessions{s,9}+5*60)));
    [~, ket_10_after] = min(abs(bins-(wt_ket_sessions{s,9}+10*60)));

    % get trans prob subsetted by type and time
%     peri_ctrl_EE = [peri_ctrl_EE; EE_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
%     peri_ctrl_EI = [peri_ctrl_EI; EI_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
%     peri_ctrl_EU = [peri_ctrl_EU; EU_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
%     peri_ctrl_IE = [peri_ctrl_IE; IE_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
%     peri_ctrl_II = [peri_ctrl_II; II_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
%     peri_ctrl_IU = [peri_ctrl_IU; IU_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
%     peri_ctrl_fromE = [peri_ctrl_fromE; fromE_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
%     peri_ctrl_fromI = [peri_ctrl_fromI; fromI_trans_prob(ctrl_5_before:ctrl_10_after,:)'];

%     peri_ket_EE = [peri_ket_EE; EE_trans_prob(ket_5_before:ket_10_after,:)'];
%     peri_ket_EI = [peri_ket_EI; EI_trans_prob(ket_5_before:ket_10_after,:)'];
%     peri_ket_EU = [peri_ket_EU; EU_trans_prob(ket_5_before:ket_10_after,:)'];
%     peri_ket_IE = [peri_ket_IE; IE_trans_prob(ket_5_before:ket_10_after,:)'];
%     peri_ket_II = [peri_ket_II; II_trans_prob(ket_5_before:ket_10_after,:)'];
%     peri_ket_IU = [peri_ket_IU; IU_trans_prob(ket_5_before:ket_10_after,:)'];
%     peri_ket_fromE = [peri_ket_fromE; fromE_trans_prob(ket_5_before:ket_10_after,:)'];
%     peri_ket_fromI = [peri_ket_fromI; fromI_trans_prob(ket_5_before:ket_10_after,:)'];

%     delta_ctrl_EE = [delta_ctrl_EE; [nanmean(EE_trans_prob(ctrl_inj:ctrl_5_after,:)-...
%         nanmean(EE_trans_prob(ctrl_5_before:ctrl_inj,:)))./nanmean(EE_trans_prob(base_start:base_end,:))]'];
%     delta_ket_EE = [delta_ket_EE; [nanmean(EE_trans_prob(ket_inj:ket_5_after,:)-...
%         nanmean(EE_trans_prob(ket_5_before:ket_inj,:)))./nanmean(EE_trans_prob(base_start:base_end,:))]'];
%     delta_ctrl_EI = [delta_ctrl_EI; [nanmean(EI_trans_prob(ctrl_inj:ctrl_5_after,:)-...
%         nanmean(EI_trans_prob(ctrl_5_before:ctrl_inj,:)))./nanmean(EI_trans_prob(base_start:base_end,:))]'];
%     delta_ket_EI = [delta_ket_EI; [nanmean(EI_trans_prob(ket_inj:ket_5_after,:)-...
%         nanmean(EI_trans_prob(ket_5_before:ket_inj,:)))./nanmean(EI_trans_prob(base_start:base_end,:))]'];
%     delta_ctrl_EU = [delta_ctrl_EU; [nanmean(EU_trans_prob(ctrl_inj:ctrl_5_after,:)-...
%         nanmean(EU_trans_prob(ctrl_5_before:ctrl_inj,:)))./nanmean(EU_trans_prob(base_start:base_end,:))]'];
%     delta_ket_EU = [delta_ket_EU; [nanmean(EU_trans_prob(ket_inj:ket_5_after,:)-...
%         nanmean(EU_trans_prob(ket_5_before:ket_inj,:)))./nanmean(EU_trans_prob(base_start:base_end,:))]'];
%     delta_ctrl_IE = [delta_ctrl_IE; [nanmean(IE_trans_prob(ctrl_inj:ctrl_5_after,:)-...
%         nanmean(IE_trans_prob(ctrl_5_before:ctrl_inj,:)))./nanmean(IE_trans_prob(base_start:base_end,:))]'];
%     delta_ket_IE = [delta_ket_IE; [nanmean(IE_trans_prob(ket_inj:ket_5_after,:)-...
%         nanmean(IE_trans_prob(ket_5_before:ket_inj,:)))./nanmean(IE_trans_prob(base_start:base_end,:))]'];
%     delta_ctrl_II = [delta_ctrl_II; [nanmean(II_trans_prob(ctrl_inj:ctrl_5_after,:)-...
%         nanmean(II_trans_prob(ctrl_5_before:ctrl_inj,:)))./nanmean(II_trans_prob(base_start:base_end,:))]'];
%     delta_ket_II = [delta_ket_II; [nanmean(II_trans_prob(ket_inj:ket_5_after,:)-...
%         nanmean(II_trans_prob(ket_5_before:ket_inj,:)))./nanmean(II_trans_prob(base_start:base_end,:))]'];
%     delta_ctrl_IU = [delta_ctrl_IU; [nanmean(IU_trans_prob(ctrl_inj:ctrl_5_after,:)-...
%         nanmean(IU_trans_prob(ctrl_5_before:ctrl_inj,:)))./nanmean(IU_trans_prob(base_start:base_end,:))]'];
%     delta_ket_IU = [delta_ket_IU; [nanmean(IU_trans_prob(ket_inj:ket_5_after,:)-...
%         nanmean(IU_trans_prob(ket_5_before:ket_inj,:)))./nanmean(IU_trans_prob(base_start:base_end,:))]'];

%     delta_ctrl_fromE = [delta_ctrl_fromE; [nanmean(fromE_trans_prob(ctrl_inj:ctrl_5_after,:)-...
%         fromE_trans_prob(ctrl_5_before:ctrl_inj,:))]'];
%     delta_ket_fromE = [delta_ket_fromE; [nanmean(fromE_trans_prob(ket_inj:ket_5_after,:)-...
%         fromE_trans_prob(ket_5_before:ket_inj,:))]'];
%     delta_ctrl_fromI = [delta_ctrl_fromI; [nanmean(fromI_trans_prob(ctrl_inj:ctrl_5_after,:)-...
%         fromI_trans_prob(ctrl_5_before:ctrl_inj,:))]'];
%     delta_ket_fromI = [delta_ket_fromI; [nanmean(fromI_trans_prob(ket_inj:ket_5_after,:)-...
%         fromI_trans_prob(ket_5_before:ket_inj,:))]'];
% 
%     peri_ctrl_fromE = [peri_ctrl_fromE; fromE_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
%     peri_ctrl_fromI = [peri_ctrl_fromI; fromI_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
%     peri_ket_fromE = [peri_ket_fromE; fromE_trans_prob(ket_5_before:ket_10_after,:)'];
%     peri_ket_fromI = [peri_ket_fromI; fromI_trans_prob(ket_5_before:ket_10_after,:)'];

    delta_ctrl_fromE = [delta_ctrl_fromE; [nanmean(fromE_trans_prob(ctrl_inj:ctrl_5_after,:)-...
        fromE_trans_prob(ctrl_5_before:ctrl_inj,:))]'];
    delta_ket_fromE = [delta_ket_fromE; [nanmean(fromE_trans_prob(ket_inj:ket_5_after,:)-...
        fromE_trans_prob(ket_5_before:ket_inj,:))]'];
    delta_ctrl_fromI = [delta_ctrl_fromI; [nanmean(fromI_trans_prob(ctrl_inj:ctrl_5_after,:)-...
        fromI_trans_prob(ctrl_5_before:ctrl_inj,:))]'];
    delta_ket_fromI = [delta_ket_fromI; [nanmean(fromI_trans_prob(ket_inj:ket_5_after,:)-...
        fromI_trans_prob(ket_5_before:ket_inj,:))]'];

    peri_ctrl_fromE = [peri_ctrl_fromE; fromE_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
    peri_ctrl_fromI = [peri_ctrl_fromI; fromI_trans_prob(ctrl_5_before:ctrl_10_after,:)'];
    peri_ket_fromE = [peri_ket_fromE; fromE_trans_prob(ket_5_before:ket_10_after,:)'];
    peri_ket_fromI = [peri_ket_fromI; fromI_trans_prob(ket_5_before:ket_10_after,:)'];

%     base_trace_fromE(count_e+1:count_e+length(fromE_connect),1:length(base_start:base_end)) = ...
%         fromE_trans_prob(base_start:base_end,:)';
%     ctrl_trace_fromE(count_e+1:count_e+length(fromE_connect),1:length(ctrl_start:ctrl_end)) = ...
%         fromE_trans_prob(ctrl_start:ctrl_end,:)';
%     ket_trace_fromE(count_e+1:count_e+length(fromE_connect),1:length(ket_start:ket_end)) = ...
%         fromE_trans_prob(ket_start:ket_end,:)';
%     base_trace_fromI(count_i+1:count_i+length(fromI_connect),1:length(base_start:base_end)) = ...
%         fromI_trans_prob(base_start:base_end,:)';
%     ctrl_trace_fromI(count_i+1:count_i+length(fromI_connect),1:length(ctrl_start:ctrl_end)) = ...
%         fromI_trans_prob(ctrl_start:ctrl_end,:)';
%     ket_trace_fromI(count_i+1:count_i+length(fromI_connect),1:length(ket_start:ket_end)) = ...
%         fromI_trans_prob(ket_start:ket_end,:)';
%     count_e = count_e + length(fromE_connect);
%     count_i = count_i + length(fromI_connect);
end

%save('\\oak-smb-giocomo.stanford.edu\groups\giocomo\fkmasuda\fkm_analysis\EAJ_revisions\transmission_probabilities.mat')

%% Plot
% 3E/F
% x = -5:10/60:10;
% plotBinnedDataOverTime(x, peri_ctrl_fromE, peri_ket_fromE)
% plotBinnedDataOverTime(x, peri_ctrl_fromI, peri_ket_fromI)
% plotBinnedDataOverTime(x, peri_ctrl_EE, peri_ket_EE)
% plotBinnedDataOverTime(x, peri_ctrl_EI, peri_ket_EI)
% plotBinnedDataOverTime(x, peri_ctrl_EU, peri_ket_EU)
% plotBinnedDataOverTime(x, peri_ctrl_IE, peri_ket_IE)
% plotBinnedDataOverTime(x, peri_ctrl_II, peri_ket_II)
% plotBinnedDataOverTime(x, peri_ctrl_IU, peri_ket_IU)

x = -4:10;
rescale = 6;
plotBinnedDataOverTime(x, bin_data(peri_ctrl_fromE,rescale), bin_data(peri_ket_fromE,rescale))
plotBinnedDataOverTime(x, bin_data(peri_ctrl_fromI,rescale), bin_data(peri_ket_fromI,rescale))

%% 3G
figure; clf; clear g;
delta_ctrl_fromE(~isfinite(delta_ctrl_fromE)) = NaN;
delta_ket_fromE(~isfinite(delta_ket_fromE)) = NaN;
delta_ctrl_fromI(~isfinite(delta_ctrl_fromI)) = NaN;
delta_ket_fromI(~isfinite(delta_ket_fromI)) = NaN;
y = vertcat(delta_ctrl_fromE,delta_ket_fromE,delta_ctrl_fromI,delta_ket_fromI);
color = vertcat(repmat({'1Ctrl_FromE'},size(delta_ctrl_fromE,1),1),...
    repmat({'2Ket_FromE'},size(delta_ket_fromE,1),1),...
    repmat({'3Ctrl_FromI'},size(delta_ctrl_fromI,1),1),...
    repmat({'4Ket_FromI'},size(delta_ket_fromI,1),1));
g=gramm('x',color,'y',y,'color',color);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',0.15);
g.draw;

%     delta_ctrl_EE,delta_ket_EE,delta_ctrl_EI,delta_ket_EI,delta_ctrl_EU,delta_ket_EU,...
%     delta_ctrl_IE,delta_ket_IE,delta_ctrl_II,delta_ket_II,delta_ctrl_IU,delta_ket_IU);
% color = vertcat(repmat({'1Ctrl_FromE'},size(delta_ctrl_fromE1,1),1),...
%     repmat({'2Ket_FromE'},size(delta_ket_fromE1,1),1),...
%     repmat({'3Ctrl_FromI'},size(delta_ctrl_fromI1,1),1),...
%     repmat({'4Ket_FromI'},size(delta_ket_fromI1,1),1),...
%     repmat({'5Ctrl_EE'},size(delta_ctrl_EE,1),1),...
%     repmat({'6Ket_EE'},size(delta_ket_EE,1),1),...
%     repmat({'7Ctrl_EI'},size(delta_ctrl_EI,1),1),...
%     repmat({'8Ket_EI'},size(delta_ket_EI,1),1),...
%     repmat({'9Ctrl_EU'},size(delta_ctrl_EU,1),1),...
%     repmat({'10Ket_EU'},size(delta_ket_EU,1),1),...
%     repmat({'11Ctrl_IE'},size(delta_ctrl_IE,1),1),...
%     repmat({'12Ket_IE'},size(delta_ket_IE,1),1),...
%     repmat({'13Ctrl_II'},size(delta_ctrl_II,1),1),...
%     repmat({'14Ket_II'},size(delta_ket_II,1),1),...
%     repmat({'15Ctrl_IU'},size(delta_ctrl_IU,1),1),...
%     repmat({'16Ket_IU'},size(delta_ket_IU,1),1));
% g=gramm('x',color,'y',y,'color',color);
% g.stat_violin('normalization','width','dodge',0,'fill','edge');
% g.stat_boxplot('width',3);
% g.set_names('x',[],'y','Transmission Probability','size',20); 
% g.set_title(sprintf('Stability Score'),'fontSize',20);
% g.axe_property('FontSize',12);
% customColorMap = [0.5 0.5 0.5g.d
%     0.8 0.2 0.8
%     0 0.8 0.2];
% g.set_color_options('map',customColorMap);
% g.draw;

%% Statistics
stats = compare_centers(delta_ctrl_fromE, delta_ket_fromE, 'paired', true);
fprintf("E: %s\n",stats.resultsstr)
stats = compare_centers(delta_ctrl_fromI, delta_ket_fromI, 'paired', true);
fprintf("I: %s\n",stats.resultsstr)