function [spdFltrTrial,spdFltrPosX,spdFltrPosT] = speedFilterData(trial,posx, post, speed, n)
% remove trial, position & time frames that have speed <n cm/s
if ~isempty(trial)
    spdFltrTrial = trial(speed>n);
else
    spdFltrTrial = [];
end

if ~isempty(posx)
    spdFltrPosX = posx(speed>n);
else
    spdFltrPosX = [];
end

if ~isempty(post)
    spdFltrPosT = post(speed>n);
else
    spdFltrPosT = [];
end


end
