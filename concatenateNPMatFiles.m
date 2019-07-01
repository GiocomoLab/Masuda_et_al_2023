function [lickt,lickx,post,posx,trial,sp] = concatenateNPMatFiles(matA, matB)
% where the data in matA occurs before matB

load(matA,'lickt','lickx','post','posx','sp','trial');
licktA = lickt; lickxA = lickx; postA = post; posxA=posx; spA = sp; trialA = trial;

load(matB,'lickt','lickx','post','posx','sp','trial');
lickt = vertcat(licktA,lickt);
lickx = vertcat(lickxA,lickx);
post = vertcat(postA,post);
posx= vertcat(posxA,posx); 
trial= vertcat(trialA,trial); 
sp = concatenateSpStructs(spA, sp);

end