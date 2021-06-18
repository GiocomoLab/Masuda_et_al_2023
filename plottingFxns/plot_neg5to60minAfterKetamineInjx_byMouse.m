function plot_neg5to60minAfterKetamineInjx_byMouse(cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Firing Rate over Time 5min before and 60 min after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleRate = 50; %hz
secInMin = 60; 
scaling  = sampleRate * secInMin / 50; %divide by 50 now that i'm downsampling in poolAllCells

figure(); clf;
hold on;

metadataCategory = 2; 
seshes = unique(cellfun(@num2str,cells.metadata(:,metadataCategory),'uni',0));
for i = 1:numel(seshes)
    
    seshIndx = ismember(cells.metadata(:,metadataCategory),seshes{i});
    xsize = size(cells.timeFRneg5to60minAfterKetamineInjx(seshIndx,:,:),2);
%     numCell = size(cells.timeFRneg5to60minAfterKetamineInjx(seshIndx,:,:),1);
    dsfactor = 1;
    x = downsample(1:xsize,dsfactor);
    x = x/scaling - 5;
    y = cells.timeFRneg5to60minAfterKetamineInjx(seshIndx,:,:);
    y = downsample(y',dsfactor)';
    y = smoothdata(y,2,'movmean',200);
    
    baselineFR = nanmean(y(:,x<0),2);
    normY = y-baselineFR;
    
    data = smoothdata(nanmean(normY,1),'gaussian',100);
    plot(x,data,'-k','LineWidth',1); hold on;
    
    postKetData = data(x>0);
    [maxValue,maxIdx] = max(postKetData);
    [minValue,minIdx] = min(postKetData);
    plot(x(maxIdx:minIdx)+find(x==0),data(maxIdx:minIdx),'-b','LineWidth',1)
    
    ketLengthTime(i) = x(minIdx+find(x==0));
    ketPeakTime(i) = x(maxIdx+find(x==0));
end



xsize = size(cells.timeFRneg5to60minAfterKetamineInjx,2);
% numCell = size(cells.timeFRneg5to60minAfterKetamineInjx,1);
dsfactor = 1;
x = downsample(1:xsize,dsfactor);
x = x/scaling - 5;
y = cells.timeFRneg5to60minAfterKetamineInjx;
y = downsample(y',dsfactor)';
y = smoothdata(y,2,'movmean',200);

baselineFR = nanmean(y(:,x<0),2);
normY = y-baselineFR;

% 
data_SEM = nanstd(normY,1)./sqrt(size(normY,1));
shadedErrorBar(x,nanmean(normY,1),data_SEM,'lineprops','-k')
plot(x,smoothdata(nanmean(normY,1),'gaussian',100),'-r','LineWidth',4)

set(gca,'TickDir','out');
set(gca,'ticklength',[0.015 0.025]);
set(gca,'layer','bottom');
% box off;
set(gca,'FontName','Helvetica');
% axis square;
set(gcf,'Position',[100 100 1000 1000])
set(gca,'FontSize',50);
% title(titleStr)
xlabel('Minutes since Injection')
ylabel('Change FR (Hz)')
xlim([-5 60])
vline(0,'-g')


end

