function [spdFltrTrial,spdFltrPosX,spdFltrPosT] = speedFilterData(trial,posx, post, speed)

% remove trial, position & time frames that have speed <n cm/s
n = 2;

spdFltrTrial = trial;
spdFltrTrial(speed<n) = [];

spdFltrPosX = posx;
spdFltrPosX(speed<n) = [];

spdFltrPosT = post;
spdFltrPosT(speed<n) = [];


end
