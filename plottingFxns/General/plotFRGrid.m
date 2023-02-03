function plotFRGrid(allCells,num2plot)
    numCells = size(allCells.spatialFR2,1);
    num2plot = min(num2plot,numCells);
    rows = ceil(sqrt(num2plot));
    cols = rows;
    figure();
    for i = 1:num2plot
 
        name = allCells.metadata{i,2};
        genotype = allCells.metadata{i,4};  
        frMap = squeeze(allCells.spatialFR2(i,:,:));

        
        subplot(rows,cols,i);
        imagesc(frMap);
        set(gca,'TickDir','out');
        set(gca,'ticklength',[0.005 0.025]);
        set(gca,'layer','bottom');
        box off;
        set(gca,'FontSize',20);
        set(gca,'FontName','Helvetica');
        set(gca,'Xticklabel',[])
        set(gca,'Yticklabel',[])
        % set(gcf,'Position',[100 100 1000 1000])
%         title(sprintf('Cell %d: %s,%s',i,name,genotype))
%         xlabel('VR cm')
%         ylabel('Trial Number')
        
    end
end