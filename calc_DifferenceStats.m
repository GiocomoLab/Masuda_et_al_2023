function [diff_means,diff_sems,p]=calc_DifferenceStats(x1,x2)

mean1 = nanmean(x1);
sem1 = std(x1)/sqrt(length(x1));

mean2 = nanmean(x2);
sem2 = std(x2)/sqrt(length(x2));

diff_means = mean1 - mean2;
diff_sems = sem1 - sem2;

[~,p] = ttest2(x1,x2);
fprintf('difference in running speed, %f ± %f cm/s; p = %f\n',diff_means,diff_sems,p)

end