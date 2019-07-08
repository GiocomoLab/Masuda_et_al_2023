sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/*.mat');
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/0Code/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));

fr_mean_allCells= nan(numel(sessions),300);
corrmatrix_allCells = nan(numel(sessions), 300, 300);

for n = 1:numel(sessions)
    try
       load(fullfile(sessions(n).folder, sessions(n).name));
       fr_mean_allCells(n,:) = mean(avg_all_fr,2);
       corrmatrix_allCells(n,:,:) = avg_all_corrmatrix;
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end

%sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/HCN1*.mat');

close all;
clear g;
x = 1:300;
y = fr_mean_allCells;

g(1,1) = gramm('x',x,'y',y);
g(1,1).stat_summary('type','std');
g(1,1).axe_property('YLim',[-0.01 0.15]);
g(1,1).set_title("Ketamine's Effect on Firing Rate (HCN1ko)", 'FontSize', 40);
g(1,1).set_names('x','Trial Number','y','Firing Rate');
g.draw()