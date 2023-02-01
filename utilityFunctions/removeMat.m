
function unitySessions = removeMat(seshID)  
    % remove .mat from session Id names
    unitySessions = extractBefore(seshID,'.mat');
end