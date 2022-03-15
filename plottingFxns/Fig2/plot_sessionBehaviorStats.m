function plot_sessionBehaviorStats(cells)
% Inputs: cells struct
% Plot Speed Behavior by sessions
% Plots speed over trials then speed over first 50 trials after each
% condition

seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));

maxTrial = 300;
session_number = numel(seshes);
avg_speed_by_trial_by_session=zeros(300,session_number);
baseline_avg_speed_by_position_by_session=zeros(session_number,400);
control_avg_speed_by_position_by_session=zeros(session_number,400);
ketamine_avg_speed_by_position_by_session=zeros(session_number,400);

for i = 1:session_number
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(cells,seshIndx);

        
    speed = extractSessionValueFromCellsStruct(seshCells.speed);
    trial = extractSessionValueFromCellsStruct(seshCells.trial);
%     post = extractSessionValueFromCellsStruct(seshCells.posT);
    posx = extractSessionValueFromCellsStruct(seshCells.posX);
%     lickt = extractSessionValueFromCellsStruct(seshCells.lickT);
%     lickx = extractSessionValueFromCellsStruct(seshCells.lickX);
    
    % Calculate the average speeds by trial for every sessions
    avg_speed_by_trial = zeros(1,maxTrial);
    for j = 1:maxTrial
        avg_speed_by_trial(j) = nanmean(speed(trial==j));
    end
    avg_speed_by_trial_by_session(:,i) = avg_speed_by_trial;
    
    % Calculate the average speeds by position for every sessions for just
    % the baseline trials
    baseline_speed = speed(trial<51);
    baseline_posx = round(posx(trial<51));
    baseline_avg_speed_by_position = zeros(1,400);
    for j = 1:400
        baseline_avg_speed_by_position(j) = nanmean(baseline_speed(baseline_posx==j));
    end
    baseline_avg_speed_by_position_by_session(i,:) = baseline_avg_speed_by_position;
    
    % Calculate the average speeds by position for every sessions for just
    % the control trials
    control_speed = speed(trial>50 & trial<101);
    control_posx = round(posx(trial>50 & trial<101));
    control_avg_speed_by_position = zeros(1,400);
    for j = 1:400
        control_avg_speed_by_position(j) = nanmean(control_speed(control_posx==j));
    end
    control_avg_speed_by_position_by_session(i,:) = control_avg_speed_by_position;
    
    % Calculate the average speeds by position for every sessions for just
    % the ketamine trials
    ketamine_speed = speed(trial>100 & trial<150);
    ketamine_posx = round(posx(trial>100 & trial<150));
    ketamine_avg_speed_by_position = zeros(1,400);
    for j = 1:400
        ketamine_avg_speed_by_position(j) = nanmean(ketamine_speed(ketamine_posx==j));
    end
    ketamine_avg_speed_by_position_by_session(i,:) = ketamine_avg_speed_by_position;
    
    


end
%% Calculate Baseline Statistics
mean_baseline_avg_speed_by_position_by_session_1 = nanmean(baseline_avg_speed_by_position_by_session(:,100:200),2);
mean_baseline_avg_speed_by_position_by_session_2 = nanmean(baseline_avg_speed_by_position_by_session(:,390:400),2);
calc_DifferenceStats(mean_baseline_avg_speed_by_position_by_session_1,mean_baseline_avg_speed_by_position_by_session_2);

%% Calculate Baseline Statistics
mean_control_avg_speed_by_position_by_session_1 = nanmean(control_avg_speed_by_position_by_session(:,100:200),2);
mean_control_avg_speed_by_position_by_session_2 = nanmean(control_avg_speed_by_position_by_session(:,390:400),2);
calc_DifferenceStats(mean_control_avg_speed_by_position_by_session_1,mean_control_avg_speed_by_position_by_session_2);

%% Calculate Baseline Statistics
mean_ketamine_avg_speed_by_position_by_session_1 = nanmean(ketamine_avg_speed_by_position_by_session(:,100:200),2);
mean_ketamine_avg_speed_by_position_by_session_2 = nanmean(ketamine_avg_speed_by_position_by_session(:,390:400),2);
calc_DifferenceStats(mean_ketamine_avg_speed_by_position_by_session_1,mean_ketamine_avg_speed_by_position_by_session_2);
%% Plot Data
close all;
clear g;
g(1,1) = gramm('x',1:290,'y',avg_speed_by_trial_by_session(1:290,:)');
g(1,1).stat_summary('setylim','true');
% g(1,1).set_title('Average Speed by Trial');
% g(1,1).set_names('x','Trials','y','Speed(cm/s)');
g(1,1).set_names('x','','y','');
g(1,1).set_color_options('map',[0 0 0]); %black

g(1,2) = gramm('x',1:400,'y',baseline_avg_speed_by_position_by_session);
g(1,2).stat_summary('setylim','true');
% g(1,2).set_title('Average Speed by Position: Baseline (50 trials)');
% g(1,2).set_names('x','Position(cm)','y','Speed(cm/s)');
g(1,2).set_names('x','','y','');
g(1,2).set_color_options('map',[0.5 0.5 0.5]); %grey
g(1,2).axe_property('YLim',[0 45]);

g(1,3) = gramm('x',1:400,'y',control_avg_speed_by_position_by_session);
g(1,3).stat_summary('setylim','true');
% g(1,3).set_title('Average Speed by Position: Post-Control(50 trials)');
% g(1,3).set_names('x','Position(cm)','y','Speed(cm/s)');
g(1,3).set_names('x','','y','');
g(1,3).set_color_options('map',[ 0.8 0.2 0.8 ]); %magenta
g(1,3).axe_property('YLim',[0 45]);

g(1,4) = gramm('x',1:400,'y',ketamine_avg_speed_by_position_by_session);
g(1,4).stat_summary('setylim','true');
% g(1,4).set_title('Average Speed by Position: Post-Ketamine(50 trials)');
% g(1,4).set_names('x','Position(cm)','y','Speed(cm/s)');
g(1,4).set_names('x','','y','');
g(1,4).set_color_options('map',[0 0.8 0.2]); %green
g(1,4).axe_property('YLim',[0 45]);
g.draw();

%%