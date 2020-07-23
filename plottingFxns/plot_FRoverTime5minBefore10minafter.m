function plot_FRoverTime5minBefore10minafter(allCells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Firing Rate over Time 5 min before injection and 10 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleRate = 50; %hz
secInMin = 60; 
scaling  = sampleRate * secInMin; 

figure(); clf;
hold on;

xsize = size(allCells.timeFRcircaControlInjx,2);
numCell = size(allCells.timeFRcircaControlInjx,1);
dsfactor = 50;
x = downsample(1:xsize,dsfactor);
x = x/scaling - 5;
y = allCells.timeFRcircaControlInjx;
plotsemline(x,y,'-m');
y = allCells.timeFRcircaKetamineInjx;
plotsemline(x,y,'-g');

set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
box off;
set(gca,'FontName','Helvetica');
axis square;
set(gcf,'Position',[100 100 1000 1000])
set(gca,'FontSize',30);
title('Control vs Ketamine')
xlabel('Minutes since Injection')
ylabel('Firing Rate (Hz)')

function plotsemline(x,y,color)
    y = downsample(y',dsfactor)';
    y = smoothdata(y,2,'movmean',30);
    data_SEM = nanstd(y,1)./sqrt(size(y,1));
    shadedErrorBar(x,nanmean(y,1),data_SEM,'lineprops',color)
    plot(x,nanmean(y,1),'-k')
end

end

