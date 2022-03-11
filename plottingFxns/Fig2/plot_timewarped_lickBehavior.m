function plot_timewarped_lickBehavior(cells)
% Plots Timewarped Lick Accuracy: 30 min post-start; 30 min post-control;
% 120min post-ketamine

addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));

seshes = unique(cellfun(@num2str,cells.metadata(:,1),'uni',0));
allLickAccuracy = nan(numel(seshes),300);


bufferMin = 10;

min = 30;
fs_scaling = 60/0.2;
time_indices = min * fs_scaling;
timewarped_lickAccuracy_postStart = nan(numel(seshes),time_indices);

time_indices_cntrl = (min+bufferMin) * fs_scaling;
timewarped_lickAccuracy_postControl = nan(numel(seshes),time_indices_cntrl);

min_ket = 120;
time_indices_ket = (min_ket+bufferMin) * fs_scaling;
timewarped_lickAccuracy_postKet = nan(numel(seshes),time_indices_ket);


for i = 1:numel(seshes)
    seshIndx = ismember(cells.metadata(:,1),seshes{i});
    seshCells = filterAllCellsStruct(cells,seshIndx);
    genotype = seshCells.metadata(:,4);
    genotypeIndx(i) = genotype(1);

    [lickAccuracyByTrial,timewarpedLickAccuracy,trial,post] = plot_SingleSessionBehavior(seshCells,false);


    postStart = find(trial==1,1);
    timewarped_lickAccuracy_postStart(i,:) = timewarpedLickAccuracy(postStart:postStart+time_indices-1);

    postControlInjx = find(trial==51,1)-(bufferMin*fs_scaling);
    timewarped_lickAccuracy_postControl(i,:) = timewarpedLickAccuracy(postControlInjx:postControlInjx+time_indices_cntrl-1);

    postKetamineInjx = find(trial==101,1)-(bufferMin*fs_scaling);
    timewarped_lickAccuracy_postKet(i,:) = timewarpedLickAccuracy(postKetamineInjx:postKetamineInjx+time_indices_ket-1);


    try
        if size(lickAccuracyByTrial,2) > 300
            allLickAccuracy(i,1:300) = lickAccuracyByTrial(1:300);
        else
            allLickAccuracy(i,:) = lickAccuracyByTrial(:);
        end
        fprintf('Finished sesssion: %d/%d\n',i,numel(seshes))
    catch
        fprintf('FAILED sesssion: %d/%d\n',i,numel(seshes)) 
    end
end
%%
% close all;
% plot_lineWithSEM(allLickAccuracy,[])
% % title('Lick Accuracy');
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% xlabel('Trials'); 
% ylabel('Lick Accuracy');
% 
% plot_lineWithSEM(timewarped_lickAccuracy_postKet,[])
% % title('Lick Accuracy');
% set(gca,'TickDir','out');
% set(gca,'ticklength',[0.005 0.025]);
% set(gca,'layer','bottom');
% box off;
% axis square;
% set(gca,'FontSize',30);
% set(gca,'FontName','Helvetica');
% xlabel('Trials'); 
% ylabel('Lick Accuracy');
%%

% for i = 1:numel(seshes)
%     seshIndx = ismember(cells.metadata(:,1),seshes{i});
%     seshCells = filterAllCellsStruct(cells,seshIndx);
%     genotype = seshCells.metadata(:,4);
%     genotypeIndx(i) = genotype(1);
% end
customColorMap = [ 0.5 0.5 0.5 %grey
    0.8 0.2 0.8 %magenta
    0 0.8 0.2]; %green

close all;
clear g;
g(1,1) = gramm('x',(0:time_indices-1)/fs_scaling,'y',timewarped_lickAccuracy_postStart);
g(1,1).stat_summary();
g(1,1).set_title('Lick Accuracy - 30min Post Start');
g(1,1).set_names('x','Time (min)','y','Lick Accuracy (% within 50cm of target)');

g(1,1).set_color_options('map',[ 0.5 0.5 0.5]); %grey


g(1,2) = gramm('x',(0:time_indices_cntrl-1)/fs_scaling-bufferMin,'y',timewarped_lickAccuracy_postControl);
g(1,2).stat_summary();
g(1,2).set_title('Lick Accuracy - 30min Post Control');
g(1,2).set_names('x','Time (min)','y','Lick Accuracy (% within 50cm of target)');
g(1,2).set_color_options('map',[ 0.8 0.2 0.8 ]); %magenta

g(1,3) = gramm('x',(0:time_indices_ket-1)/fs_scaling-bufferMin,'y',timewarped_lickAccuracy_postKet);
g(1,3).stat_summary();
g(1,3).set_title('Lick Accuracy - 120min Post Ketamine');
g(1,3).set_names('x','Time (min)','y','Lick Accuracy (% within 50cm of target)');
g(1,3).set_color_options('map',[0 0.8 0.2]); %green


g.draw();