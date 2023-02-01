%% Scratch
% scratch script for trial and error purposes 


MI_vec_sp = []; % stores the MI values of only spatial cells
MI_vec_sp_idx = []; % tells you where in cells_to_plot to find each cell
spike_depth_sp = [];

for i = 1:size(spatial_cells, 2)
    
    sp_idx = find(cells_to_plot == spatial_cells(i));
    MI_vec_sp_idx = [MI_vec_sp_idx sp_idx];
    
    this_MI = MI_vec(sp_idx);
    MI_vec_sp = [MI_vec_sp this_MI];
    
    this_depth = spike_depth(sp_idx);
    spike_depth_sp = [spike_depth_sp this_depth];
    
end

% sort MI cells based on MI score
[MI_vec_sp_sorted, MI_sort_idx] = sort(MI_vec_sp, 'descend');
spatial_cells_sorted = spatial_cells(MI_sort_idx); % these are spatial cell numbers sorted by MI

MI_vec_sp_idx_sorted = MI_vec_sp_idx(MI_sort_idx);

cell_depths = spike_depth_sp(MI_sort_idx);

% make image dir if it doesn't exist
image_save_dir = strcat('/Users/johnwen/Desktop/figs');

h = figure('Position',[100 100 160 500]); hold on; % changed from 160 to 320 for repeating, also from 500 to 166 for 50 trials

for k = 1:numel(spatial_cells_sorted)
    cla;

    fprintf('cell %d (%d/%d)\n',spatial_cells_sorted(k),k,numel(spatial_cells_sorted));

    % get spike times and index into post
    spike_t = sp.st(sp.clu==spatial_cells_sorted(k));
    [~,~,spike_idx] = histcounts(spike_t,post);

    plot(posx(spike_idx),trial(spike_idx),'k.');
    ylim([0 max(trial)+1]);
    
    xlim([0 trackLength]);
    
    title(sprintf('c%d, d=%d, MI=%1.3f',spatial_cells_sorted(k),round(cell_depths(k)),...
        MI_vec_sp_sorted(k)));
        
    % MI_vec(find(cells_to_plot == spatial_cells_sorted(k)))));
    % xticks(''); yticks('');


    % save fig
    saveas(h,fullfile(image_save_dir,sprintf('%d.png',k)),'png');
    saveas(h,fullfile(image_save_dir,sprintf('%d.pdf',k)),'pdf');

end


%% Compare cells identified as spatial across different metrics

my_cells = readtable('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/F1/F1_190620_johncontrasttrack_train1_sp_cells_by_eye.xlsx');
my_cells = table2array(my_cells);

overlap = intersect(my_cells, spatial_cells); % the cells identified by eye and MI

MI_not_mine = setdiff(spatial_cells, my_cells); % cells identified by MI but not by eye
mine_not_MI = setdiff(my_cells, spatial_cells); % cells identified by eye but not by MI


% PLOT CELLS MI Identified but I didn't
h = figure('Position',[100 100 160 500]); hold on; % changed from 160 to 320 for repeating, also from 500 to 166 for 50 trials

for k = 1:numel(MI_not_mine)
    cla;

    fprintf('cell %d (%d/%d)\n',MI_not_mine(k),k,numel(MI_not_mine));

    % get spike times and index into post
    spike_t = sp.st(sp.clu==MI_not_mine(k));
    [~,~,spike_idx] = histcounts(spike_t,post);

    plot(posx(spike_idx),trial(spike_idx),'k.');
    ylim([0 max(trial)+1]);
    
    xlim([0 trackLength]);
%     
%     title(sprintf('c%d, d=%d, MI=%1.3f',spatial_cells_sorted(k),round(cell_depths(k)),...
%         MI_vec_sp_sorted(k)));
        
    % MI_vec(find(cells_to_plot == spatial_cells_sorted(k)))));
    % xticks(''); yticks('');


    % save fig
    saveas(h,fullfile(image_save_dir2,sprintf('%d.png',k)),'png');
    saveas(h,fullfile(image_save_dir2,sprintf('%d.pdf',k)),'pdf');

end



% PLOT CELLS I IDENTIFIED BUT NOT IDENTIFIED BY MI
h = figure('Position',[100 100 160 500]); hold on; % changed from 160 to 320 for repeating, also from 500 to 166 for 50 trials

for k = 1:numel(mine_not_MI)
    cla;

    fprintf('cell %d (%d/%d)\n',mine_not_MI(k),k,numel(mine_not_MI));

    % get spike times and index into post
    spike_t = sp.st(sp.clu==mine_not_MI(k));
    [~,~,spike_idx] = histcounts(spike_t,post);

    plot(posx(spike_idx),trial(spike_idx),'k.');
    ylim([0 max(trial)+1]);
    
    xlim([0 trackLength]);
%     
%     title(sprintf('c%d, d=%d, MI=%1.3f',spatial_cells_sorted(k),round(cell_depths(k)),...
%         MI_vec_sp_sorted(k)));
        
    % MI_vec(find(cells_to_plot == spatial_cells_sorted(k)))));
    % xticks(''); yticks('');


    % save fig
    saveas(h,fullfile(image_save_dir2,sprintf('%d.png',k)),'png');
    saveas(h,fullfile(image_save_dir2,sprintf('%d.pdf',k)),'pdf');

end