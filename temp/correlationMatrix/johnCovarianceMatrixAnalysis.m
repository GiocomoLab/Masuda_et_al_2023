%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
%%

all_fr = allCellsTimeFRcircaKetamineInjx(:,:);

sampleRate = 50; %hz
min = 15;
secInMin = 60;
minBeforeInjx = 5;

preSamples = sampleRate * minBeforeInjx * secInMin;
%%
for i = 1:max(seshVector)
    % pre
    figure(1)
    datpre = zscore(all_fr(seshVector==i, 1:preSamples),0,2);
    covmatpre = cov(datpre');
    imagesc(covmatpre,[-1, 1])
    title('covmatpre');
    colorbar;

    % ket
    figure(2)
    datvr = zscore(all_fr(seshVector==i, preSamples+1:end),0,2);
    covmatket = cov(datvr');
    imagesc(covmatket,[-1, 1])
    title('covmatket');
    colorbar;
    pause
end

%%

X = categorical({'pre-Ketamine'; 'post-Ketamine'});
% X = reordercats(X,{'pre-Ketamine'; 'post-Ketamine'});
% Y = [nanmean(nanmean(covmatpre, 1), 2), nanmean(nanmean(covmatket, 1), 2)];

y = [nanmean(covmatpre, 1), nanmean(covmatket, 1)]';
numOfY = size(nanmean(covmatpre, 1),2)+size(nanmean(covmatket, 1), 2);
x = cell(numOfY,1);
x(1:size(nanmean(covmatpre, 1),2)) = {'1_preKetamine'};
x(size(nanmean(covmatpre, 1),2):end) = {'2_postKetamine'};
% Y = [nanmean(nanmean(covmatpre, 1), 2), nanmean(nanmean(covmatket, 1), 2)];
% bar(X, Y)
% title('mean covariance values vs. time')
% xlabel('portion of data')
% ylabel('mean covariance')
%%
clear g

g(1,1)=gramm('x',x,'y',abs(y),'color',x);
g(1,2)=copy(g(1));

%Averages with confidence interval
g(1,1).stat_summary('geom',{'bar','black_errorbar'},'setylim','false');
g(1,1).set_title('stat_summary()');

%Boxplots
g(1,2).stat_boxplot();
g(1,2).set_title('stat_boxplot()');

%These functions can be called on arrays of gramm objects
g.set_names('x','','y','Absolute Value of Covariance(WT)');
g.set_title('Neuronal Covariance');

figure('Position',[100 100 800 550]);
g.draw();

%%
close all;
figure(1); 
for j = 1:max(seshVector)
    clf;
    all_fr_toTest = all_fr(seshVector==j,:);

    secondsPerCovMatrx = 10;

    johnConstant = (size(all_fr_toTest,2) / sampleRate) - secondsPerCovMatrx;
    all_covScore = nan(johnConstant, 1);
    % all_covs = nan(johnConstant, size(all_fr_toTest,1), size(all_fr_toTest,1));

    % each frame in this data is 0.02
    for i = 1:johnConstant

        dat = zscore(all_fr_toTest(:, sampleRate*i:sampleRate*(i+secondsPerCovMatrx)),0,2); % check to make sure timescale isn't too large
        covmat = cov(dat');
        % corrmat = corrcov(covmat);
        % imagesc(corrmat);
    %     all_covs(i, :, :) = covmat;

%         imagesc(covmat, [-1 1])
%         set(gcf,'Position',[100 100 1000 1000])
        all_covScore(i)=nanmean(nanmean(covmat,1),2);
        % sortedCov = covmat(:, sortIdx);
    %     title(sprintf('timepoint %ss to %ss', num2str(i), num2str(i+secondsPerCovMatrx)));
    %     pause(0.005)
        % w = waitforbuttonpress


    end
    %
    plot(all_covScore)
    xlabel('time (s)')
    ylabel('mean covariance')
    title('mean covariance over time')
    pause;
end
