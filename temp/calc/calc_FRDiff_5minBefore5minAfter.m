function [cntrlDiff,ketDiff] = calc_FRDiff_5minBefore5minAfter(cells)
    sampleRate = 50; %hz
    secInMin = 60; 
    scaling  = sampleRate * secInMin; 

    cntrl = cells.timeFRcircaControlInjx;
    ket = cells.timeFRcircaKetamineInjx;
    
    cntrlBefore = cntrl(:,1:scaling*5);
    ketBefore = ket(:,1:scaling*5);
    
    cntrlAfter = cntrl(:,scaling * 5+1:scaling * 10);
    ketAfter = ket(:,scaling * 5+1:scaling * 10);
    
    cntrlDiff = mean((cntrlAfter-cntrlBefore)/nanmean(cntrlBefore),2);
    ketDiff = mean((ketAfter-ketBefore)/nanmean(cntrlBefore),2);
end

