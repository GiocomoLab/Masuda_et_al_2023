% Plot VR speed
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
% some params
params = readtable('UniversalParams.xlsx');

% calc speed
speed = calcSpeed(posx,params);

close all;
clear g;
g = gramm('x',post,'y',speed);
% g.geom_line()
g.stat_smooth()
g.draw();