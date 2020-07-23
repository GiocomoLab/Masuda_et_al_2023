
function seshIDpart = returnStrippedSeshID(seshID)
% return filename after last understore and before .mat
    newStr =  strsplit(seshID,'_');
    newStr = newStr{end};
    seshIDpart = extractBefore(newStr,'.mat');
end
