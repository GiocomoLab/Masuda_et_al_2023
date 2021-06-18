function gainModulationValues = findGainModulatedCells(allCells)
%%
nCells = size(allCells.spatialFR2,1);
gainStability = nan(nCells,2);
for j = 1:nCells
    frMap = squeeze(allCells.spatialFR2(j,:,:));


    zscore_FR = normalize(frMap,2,'center');
%         imagesc(zscore_FR)
    numTrials = size(frMap,1); 

    if numTrials >= 300
        fr1 = smoothdata(mean(frMap(281:290,:),1),'gaussian'); %pre-gain change avg
        fr2 = smoothdata(mean(frMap(291:300,:),1),'gaussian'); % post-gain change avg
        rho = corr(fr1',fr2');

%         clf;
%         figure(1); hold on;
%         plot(fr1); plot(fr2);
%         pause

    end

    gainStability(j,1) = rho;
    gainStability(j,2) = mean(fr1);
    
    close all
    plotWidth = 160;
    plotHeight = 500;
    figure('Position',[100 100 plotWidth plotHeight]); hold on;
    set(gca,'TickDir','out');
    set(gca,'ticklength',[0.005 0.025]);
    set(gca,'layer','bottom');
    box off;
    axis('tight');
    set(gca,'YDir','reverse')
    imagesc(frMap);
    title(rho);
    pause
end

