function plot_number_cells_per_mouse(cells)

% Plot  Figure

% Get Mice Name Data
metadata = cells.metadata(:,2);
% metadata = vertcat(metadata,metadata);
% metadata = metadata(~nanIndx);

g=gramm('x',metadata);

g.stat_bin('geom','stacked_bar')
g.set_title('');
g.set_names('x','Mouse','y', 'Number of Recorded Units');
customColorMap = [ 0.5 0.5 0.5 %grey
    0.8 0.2 0.8 %magenta
    0 0.8 0.2]; %green
g.set_color_options('map',customColorMap);
g.draw()
goodFigPrefs()