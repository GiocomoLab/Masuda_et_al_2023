function slope = findAttractorDynamics(cells,save_figs,sf)
% Add Description Here
nCells = size(cells.spatialFRsmooth,1);

fr = cells.spatialFRsmooth(:,1:25,:);
fr = permute(fr, [3 2 1]);
linearized_cellFR = squeeze(reshape(fr,[],1,nCells));
norm_linearized_cellFR = normalize(linearized_cellFR,1);
[rho1A,pval1] = corrcoef(norm_linearized_cellFR);
significantCellPairs = pval1<0.05;

fr = cells.spatialFRsmooth(:,26:50,:);
fr = permute(fr, [3 2 1]);
linearized_cellFR = squeeze(reshape(fr,[],1,nCells));
norm_linearized_cellFR = normalize(linearized_cellFR,1);
[rho1B,pval] = corrcoef(norm_linearized_cellFR);


fr = cells.spatialFRsmooth(:,1:50,:);
fr = permute(fr, [3 2 1]);
linearized_cellFR = squeeze(reshape(fr,[],1,nCells));
norm_linearized_cellFR = normalize(linearized_cellFR,1);
[rho1,pval] = corrcoef(norm_linearized_cellFR);

fr = cells.spatialFRsmooth(:,21:100,:);
fr = permute(fr, [3 2 1]);
linearized_cellFR = squeeze(reshape(fr,[],1,nCells));
norm_linearized_cellFR = normalize(linearized_cellFR,1);
[rho2,pval] = corrcoef(norm_linearized_cellFR);


fr = cells.spatialFRsmooth(:,101:150,:);
fr = permute(fr, [3 2 1]);
linearized_cellFR = squeeze(reshape(fr,[],1,nCells));
norm_linearized_cellFR = normalize(linearized_cellFR,1);
[rho3,pval] = corrcoef(norm_linearized_cellFR);


fr = cells.spatialFRsmooth(:,200:250,:);
fr = permute(fr, [3 2 1]);
linearized_cellFR = squeeze(reshape(fr,[],1,nCells));
norm_linearized_cellFR = normalize(linearized_cellFR,1);
[rho4,pval] = corrcoef(norm_linearized_cellFR);


fr = cells.spatialFRsmooth(:,251:290,:);
fr = permute(fr, [3 2 1]);
linearized_cellFR = squeeze(reshape(fr,[],1,nCells));
norm_linearized_cellFR = normalize(linearized_cellFR,1);
[rho5,pval] = corrcoef(norm_linearized_cellFR);


fr = cells.spatialFRsmooth(:,290:300,:);
fr = permute(fr, [3 2 1]);
linearized_cellFR = squeeze(reshape(fr,[],1,nCells));
norm_linearized_cellFR = normalize(linearized_cellFR,1);
[rho6,pval] = corrcoef(norm_linearized_cellFR);

% 
% linearized_rho1A = squeeze(reshape(rho1A,[],1));
% linearized_rho1B = squeeze(reshape(rho1B,[],1));
% 
% linearized_rho1 = squeeze(reshape(rho1,[],1));
% linearized_rho2 = squeeze(reshape(rho2,[],1));
% linearized_rho3 = squeeze(reshape(rho3,[],1));
% linearized_rho4 = squeeze(reshape(rho4,[],1));
% linearized_rho5 = squeeze(reshape(rho5,[],1));
% linearized_rho6 = squeeze(reshape(rho6,[],1));


linearized_rho1A = squeeze(reshape(rho1A,[],1));
linearized_rho1B = squeeze(reshape(rho1B,[],1));

linearized_rho1 = squeeze(reshape(rho1,[],1));
linearized_rho2 = squeeze(reshape(rho2,[],1));
linearized_rho3 = squeeze(reshape(rho3,[],1));
linearized_rho4 = squeeze(reshape(rho4,[],1));
linearized_rho5 = squeeze(reshape(rho5,[],1));
linearized_rho6 = squeeze(reshape(rho6,[],1));


% customColorMap = [ 0.5 0.5 0.5 %grey
%     0.8 0.2 0.8 %magenta
%     0 0.8 0.2]; %green

% 
% figure(1)
% tiledlayout(4,4);
% nexttile
% scatter(linearized_rho1A,linearized_rho1B)
% title("baselineA vs baselineB")
% 
% nexttile
% scatter(linearized_rho1,linearized_rho2)
% title("baseline vs cntrl")
% 
% nexttile
% scatter(linearized_rho1,linearized_rho3)
% title("baseline vs acuteKet")
% 
% nexttile
% scatter(linearized_rho2,linearized_rho3)
% title("cntrl vs acuteKet")
% 
% nexttile
% scatter(linearized_rho1,linearized_rho4)
% title("cntrl vs acuteKet")
% 
% nexttile
% scatter(linearized_rho4,linearized_rho5)
% title("acuteKet vs lateKet")
% 
% nexttile
% scatter(linearized_rho5,linearized_rho6)
% title("lateKet vs gainChange")
% 
% nexttile
% scatter(linearized_rho2,linearized_rho5)
% title("lateKet vs cntrl")
% 


%%
% threshold = 0.25;
% goodPairsIndx = (linearized_rho1A>threshold | linearized_rho1A<-threshold) & (linearized_rho1B>threshold | linearized_rho1B<-threshold);
% goodPairsIndx = linearized_rho1A>threshold | linearized_rho1A<-threshold;


linearized_significantCellPairs = squeeze(reshape(significantCellPairs,[],1));
goodPairsIndx = linearized_significantCellPairs; 
if(sum(goodPairsIndx)>10)

    clf;

    h = tiledlayout(2,4);
    nexttile
    scatter(linearized_rho1A(goodPairsIndx),linearized_rho1B(goodPairsIndx))
    title("baselineA vs baselineB")
    hl = refline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(1) = B(2);


    nexttile
    scatter(linearized_rho1(goodPairsIndx),linearized_rho2(goodPairsIndx))
    title("baseline vs cntrl")
    hl = refline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(2) = B(2);

    nexttile
    scatter(linearized_rho1(goodPairsIndx),linearized_rho3(goodPairsIndx))
    title("baseline vs acuteKet")
    hl = refline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(3) = B(2);

    nexttile
    scatter(linearized_rho2(goodPairsIndx),linearized_rho3(goodPairsIndx))
    title("cntrl vs acuteKet")
    hl = refline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(4) = B(2);

    nexttile
    scatter(linearized_rho1(goodPairsIndx),linearized_rho4(goodPairsIndx))
    title("baseline vs acuteKet")
    hl = refline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(5) = B(2);

    nexttile
    scatter(linearized_rho4(goodPairsIndx),linearized_rho5(goodPairsIndx))
    title("acuteKet vs lateKet")
    hl = refline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(6) = B(2);

    nexttile
    scatter(linearized_rho5(goodPairsIndx),linearized_rho6(goodPairsIndx))
    title("lateKet vs gainChange")
    hl = refline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(7) = B(2);

    nexttile
    scatter(linearized_rho2(goodPairsIndx),linearized_rho5(goodPairsIndx))
    title("cnrl vs lateKet")
    hl = refline;
    B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
    slope(8) = B(2);

    if save_figs
        saveas(h,fullfile(sf.image_save_dir,sprintf('%s%s%s%s%d.png',sf.name,'_',sf.sessionDate,'_sesh',sf.seshNum)),'png');
    end
end
% 
% figure(2)
% plot(slope)