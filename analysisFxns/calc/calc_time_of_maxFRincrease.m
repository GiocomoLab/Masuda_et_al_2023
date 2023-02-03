function [mean_time_of_maxFRincrease, SEM_FRtime]= calc_time_of_maxFRincrease(cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find MAX Firing Rate Time in minutes aka second derivative to find max onset of ketamine action
% input: FR: 5min before and 60 min after injection
% Found the maximum 2nd derivative (approximately the onset of ketamine) that occcurs after the ketamine
% injection within 13 minutes (half-life of ketamine) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleRate = 50; %hz
secInMin = 60; 
scaling  = sampleRate * secInMin / 50; %divide by 50 now that i'm downsampling in poolAllCells



metadataCategory = 1; %calculate by session
seshes = unique(cellfun(@num2str,cells.metadata(:,metadataCategory),'uni',0));
for i = 1:numel(seshes)
    
    seshIndx = ismember(cells.metadata(:,metadataCategory),seshes{i});
    xsize = size(cells.timeFRneg5to60minAfterKetamineInjx(seshIndx,:,:),2);
    dsfactor = 1;
    x = downsample(1:xsize,dsfactor);
    x = x/scaling - 5;
    y = cells.timeFRneg5to60minAfterKetamineInjx(seshIndx,:,:);
    y = downsample(y',dsfactor)';
    y = smoothdata(y,2,'movmean',200);
    
    baselineFR = nanmean(y(:,x<0),2);
    normY = y-baselineFR;
    
    data = smoothdata(nanmean(normY,1),'gaussian',1000);

    slope = gradient(data(:)) ./ gradient(x(:));
    curvature = gradient(slope) ./ gradient(x(:));

%     figure() %first derivavative shows maximum FR change; second postive derivative shows concave up
%     plot(x,data)
%     hold on
%     plot(x,slope)
% %     legend('y','dydx')
%     plot(x,curvature)
%     legend('y','dydx','d2ydx2')

    [maxValue,maxIdx] = max(curvature(x<13));

    time_of_maxFRincrease(i) = x(maxIdx);
end


mean_time_of_maxFRincrease = mean(time_of_maxFRincrease);
SEM_FRtime = std( time_of_maxFRincrease ) / sqrt( numel(time_of_maxFRincrease) );


end

