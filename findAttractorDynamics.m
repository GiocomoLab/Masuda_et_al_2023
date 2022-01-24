function [slope,pearson_rho] = findAttractorDynamics(cells,save_figs,sf)
% A functioning attractor network should preserve correlations between
% pairs of cells. 
% E.g. So pairs of cells that are well correlated in trial 1 through 25
% should also be well correlated in trials 26 through 50

slope = nan;
pearson_rho = nan;
%%
nCells = size(cells.spatialFRsmooth,1); %get number of cells 
% fr1to25 = cells.spatialFRsmooth(:,1:25,:); % get smoothed fr of 1st 25 trials: cellNum x trials(25) x spatial bin(200)
% fr1to25 = permute(fr1to25, [3 2 1]);
% linearized_fr1to25 = squeeze(reshape(fr1to25,[],1,nCells)); % linearize so it's a matrix of spatialbin x numCells
% [rho1A,pvalue] = corrcoef(linearized_fr1to25);
% 
% 
% fr26to50 = cells.spatialFRsmooth(:,26:50,:);
% fr26to50 = permute(fr26to50, [3 2 1]);
% linearized_fr26to50 = squeeze(reshape(fr26to50,[],1,nCells));
% [rho1B,p1b] = corrcoef(linearized_fr26to50);


fr1to50 = cells.spatialFRsmooth(:,1:50,:);% get smoothed fr of 1st 25 trials: cellNum x trials(25) x spatial bin(200)
fr1to50 = permute(fr1to50, [3 2 1]);
linearized_fr1to50 = squeeze(reshape(fr1to50,[],1,nCells));% linearize so it's a matrix of spatialbin x numCells
[rho1,pvalue] = corrcoef(linearized_fr1to50);
significantCellPairs = pvalue<0.05;  % find significant correlations
triangleMask = tril(significantCellPairs); % only grab below the diagonal to avoid repeats

fr51to100 = cells.spatialFRsmooth(:,51:100,:);
fr51to100 = permute(fr51to100, [3 2 1]);
linearized_fr51to100  = squeeze(reshape(fr51to100,[],1,nCells));
[rho2,p2] = corrcoef(linearized_fr51to100);

fr101to150 = cells.spatialFRsmooth(:,101:150,:);
fr101to150 = permute(fr101to150, [3 2 1]);
linearized_fr101to150 = squeeze(reshape(fr101to150,[],1,nCells));
[rho3,p3] = corrcoef(linearized_fr101to150);

fr200to250 = cells.spatialFRsmooth(:,200:250,:);
fr200to250 = permute(fr200to250, [3 2 1]);
linearized_fr200to250 = squeeze(reshape(fr200to250,[],1,nCells));
[rho4,p4] = corrcoef(linearized_fr200to250);

fr251to290 = cells.spatialFRsmooth(:,251:290,:);
fr251to290 = permute(fr251to290, [3 2 1]);
linearized_fr251to290 = squeeze(reshape(fr251to290,[],1,nCells));
[rho5,p5] = corrcoef(linearized_fr251to290);

fr290to300 = cells.spatialFRsmooth(:,290:300,:);
fr290to300 = permute(fr290to300, [3 2 1]);
linearized_fr290to300 = squeeze(reshape(fr290to300,[],1,nCells));
[rho6,p6] = corrcoef(linearized_fr290to300);

% 
% linearized_rho1A = squeeze(reshape(rho1A,[],1));
% linearized_rho1B = squeeze(reshape(rho1B,[],1));

linearized_significantCellPairs = squeeze(reshape(triangleMask,[],1)); % linearize mask
linearized_rho1 = squeeze(reshape(rho1,[],1));
linearized_rho2 = squeeze(reshape(rho2,[],1));
linearized_rho3 = squeeze(reshape(rho3,[],1));
linearized_rho4 = squeeze(reshape(rho4,[],1));
linearized_rho5 = squeeze(reshape(rho5,[],1));
linearized_rho6 = squeeze(reshape(rho6,[],1));


%%
% threshold = 0.25;
% goodPairsIndx = (linearized_rho1A>threshold | linearized_rho1A<-threshold) & (linearized_rho1B>threshold | linearized_rho1B<-threshold);
% goodPairsIndx = linearized_rho1A>threshold | linearized_rho1A<-threshold;

% X = linearized_rho1A(goodPairsIndx);
% Y = linearized_rho1B(goodPairsIndx);
% [~,pval] = corrcoef(X,Y,'rows', 'complete');
% significantCellPairs = pval<0.05;
% linearized_significantCellPairs = squeeze(reshape(significantCellPairs,[],1));
goodPairsIndx = linearized_significantCellPairs;



clf;
if(sum(goodPairsIndx)>10)
    h1 = figure(1);
    markerSize = 25;
    markerFaceAlpha = 0.4;
    axisSize = [-0.5 0.8 -0.5 0.8];
    tiledlayout(1,4);
    set(h1,'Position',[100 100 1600 400]);
    
%     nexttile
%     tileNum = 1;
%     X = linearized_rho1A(goodPairsIndx);
%     Y = linearized_rho1B(goodPairsIndx);
%     scatter(X,Y,markerSize,'MarkerFaceColor','k','MarkerFaceAlpha',markerFaceAlpha,'MarkerEdgeColor','w','LineWidth',0.1);
%     axis([-1 1 -1 1])
%     title("baseline vs cntrl")
%     hl = lsline;
%     B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
%     slope(tileNum) = B(2);   
%     rho = corr(X,Y,'rows', 'complete');
%     title(sprintf("baseline1 vs baseline2: rho %.3f",rho));
%     pearson_rho(tileNum) = rho;

    nexttile
    tileNum = 1;
    X = linearized_rho1(goodPairsIndx);
    Y = linearized_rho2(goodPairsIndx);
    scatter(X,Y,markerSize,'MarkerFaceColor','k','MarkerFaceAlpha',markerFaceAlpha,'MarkerEdgeColor','w','LineWidth',0.1);
    axis(axisSize)
    title("baseline vs cntrl")
    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(tileNum) = B(2);   
    rho = corr(X,Y,'rows', 'complete');
    title(sprintf("baseline vs cntrl: rho %.3f",rho));
    pearson_rho(tileNum) = rho;
    
    
    
    
    nexttile
    tileNum = 2;
    X = linearized_rho1(goodPairsIndx);
    Y = linearized_rho3(goodPairsIndx);
    scatter(X,Y,markerSize,'MarkerFaceColor','k','MarkerFaceAlpha',markerFaceAlpha,'MarkerEdgeColor','w','LineWidth',0.1);
    axis(axisSize)
    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(tileNum) = B(2);
    rho = corr(X,Y,'rows', 'complete');
    title(sprintf("baseline vs acuteKet: rho %.3f",rho));
    pearson_rho(tileNum) = rho;
    
%     nexttile
%     X = linearized_rho2(goodPairsIndx);
%     Y = linearized_rho3(goodPairsIndx);
%     scatter(X,Y)
%     title("cntrl vs acuteKet")
%     hl = lsline;
%     B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
%     slope(4) = B(2);
%     pearson_rho(4) = corr(X,Y,'rows', 'complete');

    nexttile
    tileNum = 3;
    X = linearized_rho2(goodPairsIndx);
    Y = linearized_rho4(goodPairsIndx);
    rho = corr(X,Y,'rows', 'complete');
    scatter(X,Y,markerSize,'MarkerFaceColor','k','MarkerFaceAlpha',markerFaceAlpha,'MarkerEdgeColor','w','LineWidth',0.1);
    axis(axisSize)
    title(sprintf("cntrl vs lateKet: rho %.3f",rho));
    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(tileNum) = B(2);
    pearson_rho(tileNum) = rho;

%     nexttile
%     X = linearized_rho4(goodPairsIndx);
%     Y = linearized_rho5(goodPairsIndx);
%     scatter(X,Y)
%     title("acuteKet vs lateKet")
%     hl = lsline;
%     B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
%     slope(6) = B(2);
%     pearson_rho(6) = corr(X,Y,'rows', 'complete');

    nexttile
    tileNum = 4;
    X = linearized_rho5(goodPairsIndx);
    Y = linearized_rho6(goodPairsIndx);
    scatter(X,Y,markerSize,'MarkerFaceColor','k','MarkerFaceAlpha',markerFaceAlpha,'MarkerEdgeColor','w','LineWidth',0.1);
    axis(axisSize)
    title("lateKet vs gainChange")
    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);  
    slope(tileNum) = B(2);
    rho = corr(X,Y,'rows', 'complete');
    title(sprintf("lateKet vs gainChange: rho %.3f",rho));
    pearson_rho(tileNum) = rho;
    

%     nexttile
%     X = linearized_rho2(goodPairsIndx);
%     Y = linearized_rho5(goodPairsIndx);
%     scatter(X,Y)
%     title("cnrl vs lateKet")
%     hl = lsline;
%     B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
%     slope(8) = B(2);
%     pearson_rho(8) = corr(X,Y,'rows', 'complete');
    
    if save_figs
        saveas(h1,fullfile(sf.image_save_dir,sprintf('%s%s%s%s%d.png',sf.name,'_',sf.sessionDate,'_sesh',sf.seshNum)),'png');
    end
end
