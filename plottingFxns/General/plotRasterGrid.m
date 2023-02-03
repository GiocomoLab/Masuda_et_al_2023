function plotRasterGrid(allCells,num2plot)
    numCells = size(allCells.spatialFR2,1);
    num2plot = min(num2plot,numCells);
    rows = ceil(sqrt(num2plot));
    cols = rows;
    figure();
    for i = 1:num2plot
 
        name = allCells.metadata{i,2};
        genotype = allCells.metadata{i,4};  
        lickx = allCells.lickX(i).lickx;
        lickt = allCells.lickT(i).lickt;
        posx = allCells.posX(i).posx;
        post = allCells.posT(i).post;
        trial = allCells.trial(i).trial;
        
        spike_idx = allCells.spike_idx(i);
        spike_idx = spike_idx{1};

        
        subplot(rows,cols,i);
        scatter(posx(spike_idx),trial(spike_idx),'k.');
        colormap('default')
        set(gca, 'YDir','reverse')
        ylim([0 max(trial)+1]);
        xlim([0 400]);
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