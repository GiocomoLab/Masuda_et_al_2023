function [trial, speed] = drawSpeed(data_dir)

addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
% where to find data and save images
% data_dir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190619_keicontrasttrack_ketamine1_g0';
[~,name,~] = fileparts(data_dir);
session_root = strsplit(name,'track');
session_name = {strcat(session_root{1},'track_baseline+cntrlinjx+ketamine')};

% load data
fprintf('session: %s\n',session_name{1});
load(fullfile(data_dir,strcat(session_name{1},'.mat')),'posx','post','trial'); 

% Plot VR speed
addpath(genpath('/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/UniversalParams'));
% some params
params = readtable('UniversalParams.xlsx');

% calc speed
speed = calcSpeed(posx,params);

close all;
clear g;
g = gramm('x',trial,'y',speed);
% g.geom_line()
g.stat_smooth()
g.draw();

end