function plot_stabilityTableStats(allCells, stabilityTable)
%%
addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
metadata = allCells.metadata; %session_name, animalName, sessionDate, genotype, gender, ketamine_day

%% Genotype Stairs
close all; clear g;
figure();

stabilityType = 'baselineStability';

x = stabilityTable.(stabilityType);
z = categorical(metadata(:,4));
g(1,1) = gramm('x',x,'color',z);
g(1,1).set_names('x',sprintf('%s score', stabilityType),'y','Cell Number');
% g.geom_jitter();
g(1,1).stat_bin('geom','stairs','fill','transparent');
g(1,1).set_order_options('x',{'WT','KO'},'color',{'WT','KO'});

% Genotype Violin Plot
x = categorical(metadata(:,4));
y = stabilityTable.(stabilityType);
z = categorical(metadata(:,4));
g(1,2) = gramm('x',x,'y',y,'color',z);
g(1,2).set_names('x','Genotype','y',sprintf('%s score', stabilityType));
g(1,2).stat_violin('normalization','area','dodge',0,'fill','transparent');
g(1,2).stat_boxplot('width',0.55);
g(1,2).set_order_options('x',{'WT','KO'},'color',{'WT','KO'});


% Gender
x = categorical(metadata(:,5));
y = stabilityTable.(stabilityType);
z = categorical(metadata(:,4));
g(2,1) = gramm('x',x,'y',y,'color',z);
g(2,1).set_names('x','Gender','y',sprintf('%s score', sprintf('%s score', stabilityType)));
% g.geom_jitter();
g(2,1).stat_boxplot();
g(2,1).set_order_options('color',{'WT','KO'});

% Ketamine Day
x = categorical(cell2mat(metadata(:,6)));
y = stabilityTable.(stabilityType);
z = categorical(metadata(:,4));
g(2,2) = gramm('x',x,'y',y,'color',z);
g(2,2).set_names('x','Ketamine Day','y',sprintf('%s score', stabilityType));
g(2,2).set_order_options('color',{'WT','KO'});
g(2,2).stat_boxplot();

g.set_text_options('base_size',20);
g.draw();
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
%

%% WT Only
% Genotype Stairs
close all; clear g;
figure();

stabilityType = 'baselineStability';

x = stabilityTable.(stabilityType);
z = categorical(metadata(:,4));
g(1,1) = gramm('x',x,'subset',z=='WT');
g(1,1).set_names('x',sprintf('%s score', stabilityType),'y','Cell Number');
% g.geom_jitter();
g(1,1).stat_bin('geom','stairs','fill','transparent');

% Animal Violin Plot
x = categorical(metadata(:,2));
y = stabilityTable.(stabilityType);
z = categorical(metadata(:,4));
g(1,2) = gramm('x',x,'y',y,'lightness',x,'subset',z=='WT');
g(1,2).set_names('x','Animal','y',sprintf('%s score', stabilityType));
g(1,2).stat_violin('normalization','area','dodge',0,'fill','transparent');
% g(1,2).stat_boxplot('width',0.5);


% Gender
x = categorical(metadata(:,5));
y = stabilityTable.(stabilityType);
z = categorical(metadata(:,4));
g(2,1) = gramm('x',x,'y',y,'subset',z=='WT');
g(2,1).set_names('x','Gender','y',sprintf('%s score', sprintf('%s score', stabilityType)));
% g.geom_jitter();
g(2,1).stat_boxplot();

% Ketamine Day
x = categorical(cell2mat(metadata(:,6)));
y = stabilityTable.(stabilityType);
z = categorical(metadata(:,4));
g(2,2) = gramm('x',x,'y',y,'lightness',x,'subset',z=='WT' );
g(2,2).set_names('x','Ketamine Day','y',sprintf('%s score', stabilityType));
g(2,2).stat_boxplot();

g.set_text_options('base_size',20);
g.draw();
set(gca,'FontSize',30);
set(gca,'FontName','Helvetica');
set(gcf,'Position',[100 100 1000 1000])
%
end
