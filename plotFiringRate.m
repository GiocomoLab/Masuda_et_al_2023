sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/*.mat');

%%
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/0Code/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));

fr_mean_allCells= nan(numel(sessions),300);
fr_mean_allSpatialBins = nan(numel(sessions),200);
corrmatrix_allCells = nan(numel(sessions), 300, 300);
corrblock_allCells = nan(numel(sessions), 12, 12);

for n = 1:numel(sessions)
    try
       load(fullfile(sessions(n).folder, sessions(n).name));
       fr_mean_allCells(n,:) = mean(avg_all_fr,2);
       fr_mean_allSpatialBins(n,:) = mean(avg_all_fr,1);
       corrmatrix_allCells(n,:,:) = avg_all_corrmatrix;
       corrblock_allCells(n,:,:) = avg_all_corrblock;
    catch e
        warning(e.message);
        warning('FAILED: %s\n',sessions(n).name);
    end
end

%sessions = dir('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/fkm_analysis/fr_corr_matrices/HCN1*.mat');

close all;
clear g;
x = -99:200;
y = fr_mean_allCells;

g(1,1) = gramm('x',x,'y',y);
g(1,1).stat_summary('type','std');
g(1,1).axe_property('YLim',[-0.01 0.1]);
g(1,1).set_title("Ketamine's Effect on Firing Rate by Trial(HCN1ko)", 'FontSize', 40);
g(1,1).set_names('x','Trial Number','y','Firing Rate');
g(1,1).set_text_options('base_size',20);


g.draw()


figure(2);
clims = [-0.1 0.5];
imagesc(squeeze(mean(corrblock_allCells, 1, 'omitnan')),clims);
colormap(parula(100));
colorbar;
set(gca,'XTick',0:1:12);
xticklabels(xticks*25-125)
set(gca,'YTick',0:1:12);
yticklabels(yticks*25-125)

figure(3);
clims = [-.05 0.5];
imagesc(squeeze(mean(corrmatrix_allCells, 1, 'omitnan')),clims);
colormap(parula(100));
colorbar;
set(gca,'XTick',0:10:300);
xticklabels(xticks-100)
set(gca,'YTick',0:10:300);
yticklabels(yticks-100)