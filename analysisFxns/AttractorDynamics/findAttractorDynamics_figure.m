function [slope,pearson_rho] = findAttractorDynamics_figure(cells,save_figs,sf)
% A functioning attractor network should preserve correlations between
% pairs of cells. 
% E.g. So pairs of cells that are well correlated in trial 1 through 25
% should also be well correlated in trials 26 through 50

slope = nan;
pearson_rho = nan;
nCells = size(cells.spatialFRsmooth,1); %get number of cells 

%% Find Time Binned Correlations
fr_timebinned = cells.FRtime;

trialNumStart = 1;
trialNumEnd = 50;
trialNumIndx_Start = find(cells.trial(1).trial==trialNumStart,1,'first');
trialNumIndx_End = find(cells.trial(1).trial==trialNumEnd,1,'last');
fr_time_matrix = squeeze(cell2mat(struct2cell(fr_timebinned)));
fr1to50_timebinned = fr_time_matrix(trialNumIndx_Start:trialNumIndx_End,:);% get fr of 1st 50 trials so it's a matrix of numCells x timebin
fr1to50_timebinned = smoothdata(fr1to50_timebinned,1,'gaussian',10);% smoothData
[rho1t,pvalue] = corrcoef(fr1to50_timebinned);

significantCellPairs = pvalue<0.05;  % find significant correlations
triangleMask = tril(significantCellPairs); % only grab below the diagonal to avoid repeats

trialNumStart = 51;
trialNumEnd = 100;
trialNumIndx_Start = find(cells.trial(1).trial==trialNumStart,1,'first');
trialNumIndx_End = find(cells.trial(1).trial==trialNumEnd,1,'last');
fr_time_matrix = squeeze(cell2mat(struct2cell(fr_timebinned)));
fr51to100_timebinned = fr_time_matrix(trialNumIndx_Start:trialNumIndx_End,:);% get fr of 1st 50 trials so it's a matrix of numCells x timebin
fr51to100_timebinned = smoothdata(fr51to100_timebinned,1,'gaussian',10);% smoothData
[rho2t,pvalue] = corrcoef(fr51to100_timebinned);

trialNumStart = 101;
trialNumEnd = 150;
trialNumIndx_Start = find(cells.trial(1).trial==trialNumStart,1,'first');
trialNumIndx_End = find(cells.trial(1).trial==trialNumEnd,1,'last');
fr_time_matrix = squeeze(cell2mat(struct2cell(fr_timebinned)));
fr101to150_timebinned = fr_time_matrix(trialNumIndx_Start:trialNumIndx_End,:);% get fr of 1st 50 trials so it's a matrix of numCells x timebin
fr101to150_timebinned = smoothdata(fr101to150_timebinned,1,'gaussian',10);% smoothData
[rho3t,pvalue] = corrcoef(fr101to150_timebinned);

trialNumStart = 200;
trialNumEnd = 250;
trialNumIndx_Start = find(cells.trial(1).trial==trialNumStart,1,'first');
trialNumIndx_End = find(cells.trial(1).trial==trialNumEnd,1,'last');
fr_time_matrix = squeeze(cell2mat(struct2cell(fr_timebinned)));
fr200to250_timebinned = fr_time_matrix(trialNumIndx_Start:trialNumIndx_End,:);% get fr of 1st 50 trials so it's a matrix of numCells x timebin
fr200to250_timebinned = smoothdata(fr200to250_timebinned,1,'gaussian',10);% smoothData
[rho4t,pvalue] = corrcoef(fr200to250_timebinned);

trialNumStart = 251;
trialNumEnd = 290;
trialNumIndx_Start = find(cells.trial(1).trial==trialNumStart,1,'first');
trialNumIndx_End = find(cells.trial(1).trial==trialNumEnd,1,'last');
fr_time_matrix = squeeze(cell2mat(struct2cell(fr_timebinned)));
fr251to290_timebinned = fr_time_matrix(trialNumIndx_Start:trialNumIndx_End,:);% get fr of 1st 50 trials so it's a matrix of numCells x timebin
fr251to290_timebinned = smoothdata(fr251to290_timebinned,1,'gaussian',10);% smoothData
[rho5t,pvalue] = corrcoef(fr251to290_timebinned);

trialNumStart = 291;
trialNumEnd = 300;
trialNumIndx_Start = find(cells.trial(1).trial==trialNumStart,1,'first');
trialNumIndx_End = find(cells.trial(1).trial==trialNumEnd,1,'last');
fr_time_matrix = squeeze(cell2mat(struct2cell(fr_timebinned)));
fr291to300_timebinned = fr_time_matrix(trialNumIndx_Start:trialNumIndx_End,:);% get fr of 1st 50 trials so it's a matrix of numCells x timebin
fr291to300_timebinned = smoothdata(fr291to300_timebinned,1,'gaussian',10);% smoothData
[rho6t,pvalue] = corrcoef(fr291to300_timebinned);


timebinned_significantCellPairs = squeeze(reshape(triangleMask,[],1)); % linearize mask
timebinned_rho1t = squeeze(reshape(rho1t,[],1));
timebinned_rho2t = squeeze(reshape(rho2t,[],1));
timebinned_rho3t = squeeze(reshape(rho3t,[],1));
timebinned_rho4t = squeeze(reshape(rho4t,[],1));
timebinned_rho5t = squeeze(reshape(rho5t,[],1));
timebinned_rho6t = squeeze(reshape(rho6t,[],1));
%%

goodPairsIndx = timebinned_significantCellPairs;



clf;
if(sum(goodPairsIndx)>5)
    h1 = figure(1);
    markerSize = 25;
    markerFaceAlpha = 0.4;
    axisSize = [-0.5 0.8 -0.5 0.8];
    tiledlayout(1,4);

    set(h1,'Position',[100 100 1600 400]);
%%%% 
    % Comparing timebinned binned
    
    nexttile
    tileNum = 1;
    X = timebinned_rho1t(timebinned_significantCellPairs);
    Y = timebinned_rho2t(timebinned_significantCellPairs);
    scatter(X,Y,markerSize,'MarkerFaceColor','k','MarkerFaceAlpha',markerFaceAlpha,'MarkerEdgeColor','w','LineWidth',0.1);
    axis(axisSize)
    title("baseline(time) vs cntrl(time)")
    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(tileNum) = B(2);   
    rho = corr(X,Y,'rows', 'complete');
    title(sprintf("baseline vs cntrl: rho %.3f",rho));
    pearson_rho(tileNum) = rho;
    
 
    nexttile
    tileNum = 2;
    X = timebinned_rho1t(timebinned_significantCellPairs);
    Y = timebinned_rho3t(timebinned_significantCellPairs);
    scatter(X,Y,markerSize,'MarkerFaceColor','k','MarkerFaceAlpha',markerFaceAlpha,'MarkerEdgeColor','w','LineWidth',0.1);
    axis(axisSize)
    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(tileNum) = B(2);
    rho = corr(X,Y,'rows', 'complete');
    title(sprintf("baseline(time) vs acuteKet(time): rho %.3f",rho));
    pearson_rho(tileNum) = rho;
   

    nexttile
    tileNum = 3;
    X = timebinned_rho1t(timebinned_significantCellPairs);
    Y = timebinned_rho5t(timebinned_significantCellPairs);
    rho = corr(X,Y,'rows', 'complete');
    scatter(X,Y,markerSize,'MarkerFaceColor','k','MarkerFaceAlpha',markerFaceAlpha,'MarkerEdgeColor','w','LineWidth',0.1);
    axis(axisSize)
    title(sprintf("baseline(time) vs lateKet(time): rho %.3f",rho));
    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(tileNum) = B(2);
    pearson_rho(tileNum) = rho;

    nexttile
    tileNum = 4;
    X = timebinned_rho5t(timebinned_significantCellPairs);
    Y = timebinned_rho6t(timebinned_significantCellPairs);
    scatter(X,Y,markerSize,'MarkerFaceColor','k','MarkerFaceAlpha',markerFaceAlpha,'MarkerEdgeColor','w','LineWidth',0.1);
    axis(axisSize)
    hl = lsline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);  
    slope(tileNum) = B(2);
    rho = corr(X,Y,'rows', 'complete');
    title(sprintf("lateKet(time) vs gainChange(time): rho %.3f",rho));
    pearson_rho(tileNum) = rho;


    
    if save_figs
        saveas(h1,fullfile(sf.image_save_dir,sprintf('%s%s%s%s%d.png',sf.name,'_',sf.sessionDate,'_sesh',sf.seshNum)),'png');
    end
end
