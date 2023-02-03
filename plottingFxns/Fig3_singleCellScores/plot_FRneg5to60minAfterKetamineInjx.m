function plot_FRneg5to60minAfterKetamineInjx(cells,titleStr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Firing Rate over Time 5min before and 60 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleRate = 50; %hz
secInMin = 60; 
scaling  = sampleRate * secInMin / 50; %divide by 50 now that i'm downsampling in poolAllCells

figure(); clf;
hold on;

xsize = size(cells.timeFRneg5to60minAfterKetamineInjx,2);
numCell = size(cells.timeFRneg5to60minAfterKetamineInjx,1);
dsfactor = 1;
x = downsample(1:xsize,dsfactor);
x = x/scaling - 5;
y = cells.timeFRneg5to60minAfterKetamineInjx;
y = downsample(y',dsfactor)';
y = smoothdata(y,2,'movmean',100);
data_SEM = nanstd(y,1)./sqrt(size(y,1));

shadedErrorBar(x,nanmean(y,1),data_SEM,'lineprops','-k')
plot(x,nanmean(y,1),'-k','LineWidth',2)
set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
axis square;
set(gcf,'Position',[100 100 1000 1000])
set(gca,'FontSize',30);
title(titleStr)
xlabel('Minutes since Injection')
ylabel('Firing Rate (Hz)')
xlim([-5 60])


end

