function plot_lineWithSEM(matrix,indx)
    x = 1:size(matrix,2);
    if isempty(indx)
        y = nanmean(matrix,1);
    else
        y = nanmean(matrix(indx,:),1);
    end
    data_SEM = nanstd(matrix,1)./sqrt(size(matrix,1));
    shadedErrorBar(x,y,data_SEM)
    
    set(gca,'TickDir','out');
    set(gca,'ticklength',[0.015 0.025]);
    set(gca,'layer','bottom');
    box off;
    set(gca,'FontName','Helvetica');
end