function lickAccuracyByTrial = plot_SingleSessionBehavior(cells,savefigTF)

    
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
% load data
trial = extractSessionValueFromCellsStruct(cells.trial);
post = extractSessionValueFromCellsStruct(cells.posT);
posx = extractSessionValueFromCellsStruct(cells.posX);
post = extractSessionValueFromCellsStruct(cells.posT);
posx = extractSessionValueFromCellsStruct(cells.posX);
lickt = extractSessionValueFromCellsStruct(cells.lickT);
lickx = extractSessionValueFromCellsStruct(cells.lickX);

session_name = sprintf('%s_%s_%s_%s', cells.metadata{1,2},cells.metadata{1,3},cells.metadata{1,4},cells.metadata{1,8});
fprintf('session: %s\n',session_name);


[~,~,lick_idx] = histcounts(lickt,post);

lickAccuracyByTrial = zeros(1,max(trial));
for i = 1:max(trial)
    trialLicks = lickx(trial(lick_idx) == i);
    goodLicks = sum(trialLicks<25) + sum(trialLicks>max(posx)-25); 
    if trialLicks ~= 0
        lickAccuracyByTrial(i) = goodLicks/numel(trialLicks);
    else
        lickAccuracyByTrial(i) = 0.0;
    end
end


%%
if savefigTF

    clear g;
    figure('Position',[100 100 1800 600]);
    g(1,1) = gramm('y',trial(lick_idx)*-1,'x',posx(lick_idx));
    g(1,1).geom_point();
    g(1,1).set_title(session_name);
    g(1,1).set_names('x','Lick Position in VR','y','trial number');

    g(1,2) = gramm('y',trial(lick_idx),'x',posx(lick_idx));
    g(1,2).stat_bin();
    g(1,2).set_title('Binned Licks by Position');
    g(1,2).set_names('x','Lick Position in VR','y','Lick Number');

    g(1,3) = gramm('x',trial(lick_idx),'y',posx(lick_idx));
    g(1,3).stat_bin('nbins',30,'geom','line');
    g(1,3).set_title('Binned Licks by Trial');
    g(1,3).set_names('x','Trial Number','y','Number of Licks per 10 Trials');

    g(1,4) = gramm('x',1:max(trial),'y',lickAccuracyByTrial);
    g(1,4).stat_smooth();
    g(1,4).set_title('Lick Accuracy - Smoothed');
    g(1,4).set_names('x','Trial Number','y','Lick Accuracy (% within 50cm of target)');
    g.draw();




    imgDir = '/Users/KeiMasuda/Desktop/fkm_analysis/behavior';
    saveName = fullfile(imgDir, strcat(session_name,'_behavior.jpg'));
    saveas(gcf,saveName, 'jpg'); 
end

