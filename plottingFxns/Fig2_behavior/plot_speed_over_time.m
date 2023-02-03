function plot_speed_over_time(cells)
% Plots Speed over Time: 30 min post-start; 30 min post-control;
% 120min post-ketamine

addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));

seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));
% allspeed = nan(numel(seshes),300);

bufferMin = 0;
min = 30;
min_ket = 30;
fs_scaling = 60/0.2;

% preallocate vector space for each epoch
time_indices = min * fs_scaling;
timewarped_speed_postStart = nan(numel(seshes),time_indices);

time_indices_cntrl = (min+bufferMin) * fs_scaling;
timewarped_speed_postControl = nan(numel(seshes),time_indices_cntrl);

time_indices_ket = (min_ket+bufferMin) * fs_scaling;
timewarped_speed_postKet = nan(numel(seshes),time_indices_ket);

% preallocate vector space for each epoch
for i = 1:numel(seshes)
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(cells,seshIndx);

    speed = extractSessionValueFromCellsStruct(seshCells.speed);
    trial = extractSessionValueFromCellsStruct(seshCells.trial);

    postStart = find(trial==1,1);
    timewarped_speed_postStart(i,:) = speed(postStart:postStart+time_indices-1);

    postControlInjx = find(trial==51,1)-(bufferMin*fs_scaling);
    timewarped_speed_postControl(i,:) = speed(postControlInjx:postControlInjx+time_indices_cntrl-1);

    postKetamineInjx = find(trial==101,1)-(bufferMin*fs_scaling);
    timewarped_speed_postKet(i,:) = speed(postKetamineInjx:postKetamineInjx+time_indices_ket-1);
    
end
%%
% Plot Speed Graphs on Top of Each other (requires min to be all the same for all 3 epochs)
customColorMap = [ 0.5 0.5 0.5 %grey
    0.8 0.2 0.8 %magenta
    0 0.8 0.2]; %green

% Concatenate data vectors
x = (0:time_indices-1)/fs_scaling-bufferMin;
y = vertcat(timewarped_speed_postStart, timewarped_speed_postControl,timewarped_speed_postKet);
colors = horzcat(repmat({'Baseline'},1,size(timewarped_speed_postStart,1)),repmat({'Control'},1,size(timewarped_speed_postControl,1)),repmat({'Ketamine'},1,size(timewarped_speed_postKet,1)));

close all;
figure(1)
g = gramm('x',x,'y',smoothdata(y,2),'color',colors);
g.stat_summary('setylim','true');
g.set_title(sprintf('Speed (cm/s) - %d min',min));
g.set_names('x','Time (min)','y','Speed (cm/s)');
g.set_color_options('map',customColorMap); %magenta
% g.axe_property('YLim',[0 50]);
g.draw


% g.stat_summary('geom',{'edge_bar','black_errorbar'},'type','ci','setylim','true');
%%
% Plot Each speed graph separately 
close all;
figure(2)
clear g;
g(1,1) = gramm('x',(0:time_indices-1)/fs_scaling,'y',smoothdata(timewarped_speed_postStart,2));
g(1,1).stat_summary('bin_in',1000,'setylim','true');
g(1,1).set_title('Speed (cm/s) - 30min Post Start');
g(1,1).set_names('x','Time (min)','y','Speed (cm/s)');
g(1,1).set_color_options('map',[ 0.5 0.5 0.5]); %grey
% g.axe_property('YLim',[0 50]);
% g.draw

g(1,2) = gramm('x',(0:time_indices_cntrl-1)/fs_scaling-bufferMin,'y',smoothdata(timewarped_speed_postControl,2));
g(1,2).stat_summary('setylim','true');
g(1,2).set_title('Speed (cm/s) - 30min Control');
g(1,2).set_names('x','Time (min)','y','Speed (cm/s)');
g(1,2).set_color_options('map',[0.8 0.2 0.8]); %magenta
% g.axe_property('YLim',[0 50]);
% g.draw

g(1,3) = gramm('x',(0:time_indices_ket-1)/fs_scaling-bufferMin,'y',smoothdata(timewarped_speed_postKet,2));
g(1,3).stat_summary('setylim','true');
% g(1,1).geom_line();
g(1,3).set_title('Speed (cm/s) - 120min Ketamine');
g(1,3).set_names('x','Time (min)','y','Speed (cm/s)');
g(1,3).set_color_options('map',[0 0.8 0.2]); %green
% g.axe_property('YLim',[0 50]);
g.draw();

