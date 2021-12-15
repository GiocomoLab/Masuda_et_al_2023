function [cntrlDiff,ketDiff] = find_frDiff_5minBefore5minafter(cells)
sampleRate = 50; %hz
secInMin = 60; 
scaling  = sampleRate * secInMin; 
% dsfactor = 50;
%%
cntrl = cells.timeFRcircaControlInjx;
ket = cells.timeFRcircaKetamineInjx;

cntrlBefore = cntrl(:,1:scaling*5);
ketBefore = ket(:,1:scaling*5);

cntrlAfter = cntrl(:,scaling * 5+1:scaling * 10);
ketAfter = ket(:,scaling * 5+1:scaling * 10);

% cntrlDiff = mean(cntrlAfter-cntrlBefore,2);
% ketDiff = mean(ketAfter-ketBefore,2);

cntrlDiff = mean((cntrlAfter-cntrlBefore)./nanmean(cntrlBefore),2);
ketDiff = mean((ketAfter-ketBefore)./nanmean(ketBefore),2);