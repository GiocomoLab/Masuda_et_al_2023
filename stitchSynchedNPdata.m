function [lickt,lickx,post,posx,sp,trial,trial_contrast, trial_gain] =...
    stitchSynchedNPdata(data_dir, session_name, baseline_name, controlinjx_name, drug_name)

% % where to find data and save images
data_dir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190625_keicontrasttrack_propofol1_g0';
session_name = 'G4_190625_keicontrasttrack_baseline+cntrlinjx+propofol';
baseline_name  = 'G4_190625_keicontrasttrack_baseline1';
controlinjx_name = 'G4_190625_keicontrasttrack_controlinjx1';
drug_name = 'G4_190625_keicontrasttrack_propofol1';

% load data
fprintf('session: %s\n',session_name);
load(fullfile(data_dir,strcat(baseline_name,'.mat')),'lickt','lickx','post','posx','sp','trial');%load baseline mat file
licktCombo = lickt;
lickxCombo = lickx;
postCombo = post;
posxCombo=posx; 
spCombo = sp;
trialCombo = trial;

 
load(fullfile(data_dir,strcat(controlinjx_name,'.mat')),'lickt','lickx','post','posx','sp','trial');%load controlinjx mat file
concatenateSpStructs

load(fullfile(data_dir,strcat(drug_name,'.mat')),'lickt','lickx','post','posx','sp','trial');%load drug mat file




load(fullfile(data_dir,strcat(session_name1{1},'.mat')),'lickt','lickx','post','posx','sp','trial'); %presaline 0

end

function comboMatFile = concatenateNPMatFiles(matA, matB)
% where the data in structA occurs before structB

comboStruct = structA;
comboStruct.st = horzcat(structA.st, structB.st);
comboStruct.spikeTemplates = horzcat(structA.spikeTemplates, structB.spikeTemplates);
comboStruct.clu = horzcat(structA.clu, structB.clu);
comboStruct.tempScalingAmps = horzcat(structA.tempScalingAmps, structB.tempScalingAmps);


end


function comboStruct = concatenateSpStructs(structA, structB)
% where the data in structA occurs before structB

comboStruct = structA;
comboStruct.st = horzcat(structA.st, structB.st);
comboStruct.spikeTemplates = horzcat(structA.spikeTemplates, structB.spikeTemplates);
comboStruct.clu = horzcat(structA.clu, structB.clu);
comboStruct.tempScalingAmps = horzcat(structA.tempScalingAmps, structB.tempScalingAmps);


end