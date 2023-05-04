function plot_STATS_5minBefore5minafter(cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats comparing Firing Rate over Time 5 min before injection and 5 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleRate = 50; %hz
secInMin = 60; 
scaling  = sampleRate * secInMin; 
% dsfactor = 50;
%
cntrl = cells.timeFRcircaControlInjx;
ket = cells.timeFRcircaKetamineInjx;

cntrlBefore = cntrl(:,1:scaling*5);
ketBefore = ket(:,1:scaling*5);

cntrlAfter = cntrl(:,scaling * 5+1:scaling * 10);
ketAfter = ket(:,scaling * 5+1:scaling * 10);

% cntrlDiff = mean(cntrlAfter,2)-mean(cntrlBefore,2);
% ketDiff = mean(ketAfter,2)-mean(ketBefore,2);

% cntrlDiff = mean(cntrlAfter-cntrlBefore,2);
% ketDiff = mean(ketAfter-ketBefore,2);
% 
cntrlDiff = mean((cntrlAfter-cntrlBefore)/nanmean(cntrlBefore),2);
ketDiff = mean((ketAfter-ketBefore)/nanmean(cntrlBefore),2);
fprintf('\nCntrl Diff: %f',mean(mean((cntrlAfter-cntrlBefore))))
fprintf('\nKet Diff: %f',mean(mean((ketAfter-ketBefore))))

% absCntrlDiff = abs(cntrlDiff);
% absKetDiff = abs(cntrlDiff);

[h,p] = ttest2(cntrlDiff, ketDiff);
fprintf('\n')
fprintf(num2str(p))
p = ranksum(cntrlDiff,ketDiff);
fprintf('\n')
fprintf(num2str(p))

% ket_sem = std( ketDiff ) / sqrt( numel(ketDiff) );
% cntrl_sem = std( cntrlDiff ) / sqrt( numel(cntrlDiff) );
%%
clear g;
figure(); hold on;
y = vertcat(cntrlDiff,ketDiff);

x = vertcat(...
    repmat({'Control'}, size(cntrlDiff)),...
    repmat({'Ketamine'}, size(ketDiff))...
    );

boxplot(y,x,'ColorGroup',x,'colors','k','symbol','','PlotStyle','traditional','Widths',0.5)
set(gca,'TickDir','out');
set(gca,'ticklength',[0.005 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
ylabel('Hz')
set(findobj(gca,'type','line'),'linew',2)
title(sprintf('P-value: %0.9f',p))

%% plot all values over sessions
ctrl_sum = 0;
ket_sum = 0;
counter = 0;
sess_idx = 1;
curr_sess = cells.metadata{1,1};
for c = 1:length(cells.metadata)
	if ~strcmp(cells.metadata{c,1},curr_sess)
		ctrl_diff(sess_idx) = ctrl_sum/counter;
		ket_diff(sess_idx) = ket_sum/counter;
		sess_idx = sess_idx + 1;
		curr_sess = cells.metadata{c,1};
		ctrl_sum = 0;
		ket_sum = 0;
		counter = 0;
	end
	ctrl_sum = ctrl_sum + cntrlDiff(c);
	ket_sum = ket_sum + ketDiff(c);
	counter = counter + 1;
end
ctrl_diff(sess_idx+1) = ctrl_sum/counter;
ket_diff(sess_idx+1) = ket_sum/counter;

figure; clf; clear g;
y = vertcat(ctrl_diff',ket_diff');
x = vertcat(...
    repmat({'Control'}, size(ctrl_diff')),...
    repmat({'Ketamine'}, size(ket_diff'))...
    );
g=gramm('x',x,'y',y,'color',x);
g.stat_violin('normalization','width','dodge',0,'fill','edge');
g.stat_boxplot('width',0.15);
g.draw;
ylim([-1 2])
end

