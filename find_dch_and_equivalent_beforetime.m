function dchAndBefore = find_dch_and_equivalent_beforetime(cells)

dsfactor = 100;

% initialize variables
dch_avg_cellFR = nan(size(cells.FRtime));
before_dch_avg_cellFR = nan(size(cells.FRtime));

for i = 1:numel(cells.FRtime)

    fr_time = cells.FRtime(i).FRtime;
    smoothedCellFR = smoothdata(fr_time, 'gaussian',50);
    downsample_fr_time = downsample(smoothedCellFR,dsfactor);
    dch = cells.dch(i).dch;
    decoherenceTimeIdx_correctForFs = dch.decoherenceTimeIdx*dch.Fs; %decoherenceTimeIdx is literal time so need to multiply by Fs to get indx
    
    if ~isempty(decoherenceTimeIdx_correctForFs)
        dch_fr_time = downsample_fr_time(decoherenceTimeIdx_correctForFs); 
        dch_avg_cellFR(i) = mean(dch_fr_time);

        dchStartIndx = decoherenceTimeIdx_correctForFs(1);
        dchIndxLength = numel(decoherenceTimeIdx_correctForFs);
        beforeStartIndx = dchStartIndx-dchIndxLength+1;
%         beforeStartIndx = 1:dchIndxLength;
        if beforeStartIndx < 1
            beforeStartIndx = 1;
        end
        before_dch_fr_time = downsample_fr_time(beforeStartIndx:dchStartIndx);
        before_dch_avg_cellFR(i) = mean(before_dch_fr_time);
    end

end

dchAndBefore.before_dch_avg_cellFR = before_dch_avg_cellFR;
dchAndBefore.dch_avg_cellFR = dch_avg_cellFR;