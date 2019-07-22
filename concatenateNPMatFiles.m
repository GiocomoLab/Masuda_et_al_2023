function [lickt,lickx,post,posx,trial,sp] = concatenateNPMatFiles(matA, matB)
% where the data in matA occurs before matB

load(matA,'lickt','lickx','post','posx','sp','trial');
if isempty(post) || isempty(posx) || isempty(sp) || isempty(trial)
    fprintf(strcat('Coud not find data in: \n', matA, '\n'));
    return
end
licktA = lickt; lickxA = lickx; postA = post; posxA=posx; spA = sp; trialA = trial;


load(matB,'lickt','lickx','post','posx','sp','trial');
if isempty(post) || isempty(posx) || isempty(sp) || isempty(trial)
    fprintf(strcat('Coud not find data in: \n', matB, '\n'));
    return
end

lickt = vertcat(licktA,licktA(end)+lickt);
lickx = vertcat(lickxA,lickx);
post = vertcat(postA,postA(end)+post);
posx= vertcat(posxA,posx); 
trial= vertcat(trialA,trialA(end)+trial); 
sp = concatenateSpStructs(spA, sp);

end