function [maxFRtime, SEM_FRtime]= calc_maxFRtime(cells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find MAX Firing Rate  Time in minutes
% input: FR: 5min before and 60 min after injection
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
    
    data = smoothdata(nanmean(normY,1),'gaussian',100);


    [maxValue,maxIdx] = max(data);
    [minValue,minIdx] = min(data);
    
    ketLengthTime(i) = x(minIdx);
    ketPeakTime(i) = x(maxIdx);
end


maxFRtime = mean(ketPeakTime);
SEM_FRtime = std( ketPeakTime ) / sqrt( numel(ketPeakTime) );


end

