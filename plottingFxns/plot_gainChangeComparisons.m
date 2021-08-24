
% [drugFRdiff,cntrlFRdiff, drugFREffectScore, drugCorrEffectScore, spike_depth(k)];

% plottting gainModulation vs Drug Effect
scatter(wt_ket_Cells.drugEffectScores(:,3),wt_ket_Cells.gainModulationValues(:,1))
scatter(wt_ket_Cells_onlyStable.drugEffectScores(:,3),wt_ket_Cells_onlyStable.gainModulationValues(:,1))

% plottting gainModulation vs drugFRdiff
scatter(wt_ket_Cells_onlyStable.drugEffectScores(:,1),wt_ket_Cells_onlyStable.gainModulationValues(:,1))
scatter(wt_ket_Cells.drugEffectScores(:,1),wt_ket_Cells.gainModulationValues(:,1))
scatter(abs(wt_ket_Cells.drugEffectScores(:,1)),wt_ket_Cells.gainModulationValues(:,1))

% looking at drugFREffectScore
scatter(wt_ket_Cells.drugEffectScores(:,4),wt_ket_Cells.gainModulationValues(:,1))
scatter(wt_ket_Cells.drugEffectScores(:,4),abs(wt_ket_Cells.gainModulationValues(:,1)))
scatter(abs(wt_ket_Cells.drugEffectScores(:,4)),abs(wt_ket_Cells.gainModulationValues(:,1)))

scatter(wt_ket_Cells.drugEffectScores(:,2),wt_ket_Cells.gainModulationValues(:,1))
scatter(wt_ket_Cells.drugEffectScores(:,2),abs(wt_ket_Cells.gainModulationValues(:,1)))
scatter(abs(wt_ket_Cells.drugEffectScores(:,2)),abs(wt_ket_Cells.gainModulationValues(:,1)))

scatter(wt_ket_Cells.drugEffectScores(:,3),wt_ket_Cells.gainModulationValues(:,1))
scatter(wt_ket_Cells.drugEffectScores(:,3),abs(wt_ket_Cells.gainModulationValues(:,1)))
%%
gainModValues = wt_ket_Cells_onlyStable.gainModulationValues(:,1);
gainNoNans = gainModValues(~isnan(gainModValues));
drugEffectScores_NoGainNans = wt_ket_Cells_onlyStable.drugEffectScores(~isnan(gainModValues),:);
scatter(zscore(drugEffectScores_NoGainNans(:,3)),zscore(gainNoNans))

%%
metadata = allCells.metadata; %session_name, animalName, sessionDate, genotype, gender, ketamine_day


% Genotype Stairs
close all; clear g;
figure();

g = gramm('x',gainNoNans,'color',c);
g.stat_bin('geom','stairs','fill','transparent');
% g.set_order_options('x',{'WT','KO'},'color',{'WT','KO'});
g.draw(); 

% plottting gainModulation vs Drug Effect
scatter(allCells.drugEffectScores(:,3),allCells.gainModulationValues(:,1))
scatter(allCells.drugEffectScores(:,3),allCells.gainModulationValues(:,1))