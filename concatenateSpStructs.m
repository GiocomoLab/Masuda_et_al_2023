function comboStruct = concatenateSpStructs(structA, structB)
% where the data in structA occurs before structB

comboStruct = structA;
comboStruct.st = vertcat(structA.st, structB.st);
comboStruct.spikeTemplates = vertcat(structA.spikeTemplates, structB.spikeTemplates);
comboStruct.clu = vertcat(structA.clu, structB.clu);
comboStruct.tempScalingAmps = vertcat(structA.tempScalingAmps, structB.tempScalingAmps);


end