function stabilityTable = findStableCells(allCells)
%%
    allCellsFR = allCells.spatialFR;
    nCells = size(allCellsFR,1);
    stabilityValues = nan(nCells,5);
    for k = 1:nCells

    %     frMap = squeeze(allCells.spatialFR(k,:,:));
    %     close all
    %     plotWidth = 160;
    %     plotHeight = 500;
    %     figure('Position',[100 100 plotWidth plotHeight]); hold on;
    %     set(gca,'TickDir','out');
    %     set(gca,'ticklength',[0.005 0.025]);
    %     set(gca,'layer','bottom');
    %     box off;
    %     axis('tight');
    %     set(gca,'YDir','reverse')
    %     imagesc(frMap);
    %     
        stabilityScoreCurve = allCells.stabilityScoreCurve(k,:);
    %     figure()
    %     plot(stabilityScoreCurve)

        totalStability = nanmean(stabilityScoreCurve(1,:));
        baselineStability = nanmean(stabilityScoreCurve(1,1:50));
        acuteDrugStability = nanmean(stabilityScoreCurve(1,100:150));
        endingStability = nanmean(stabilityScoreCurve(1,190:290));
        gainStability = nanmean(stabilityScoreCurve(1,290:300));


        stabilityValues(k,:) = [totalStability, baselineStability, acuteDrugStability, endingStability, gainStability];
    end

    stabilityTable = array2table(stabilityValues, 'VariableNames',...
        {'totalStability', 'baselineStability', 'acuteDrugStability', 'endingStability', 'gainStability'});
    
    
    
end
