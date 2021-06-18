function cells_autocorrelationScore = calc_dch_autocorrelationScore_old(cells)
close all;
cells_autocorrelationScore = nan(numel(cells.metadata),1);
for i = 1:size(cells.metadata,1)
   
    smoothFR = cells.spatialFRsmooth;
    try
        dch = cells.dch(i).dch;
        decoherenceTrialIdx = dch.decoherenceIdx;
    end
    
    if ~isempty(decoherenceTrialIdx)
        dchIndxLength = numel(decoherenceTrialIdx);
        controlTrials = 101-dchIndxLength:100;
        if controlTrials(1) < 1
            controlTrials = 1:dchIndxLength;
        end
        control_spatialFR = squeeze(smoothFR(i,controlTrials,:));
        dch_spatialFR = squeeze(smoothFR(i,decoherenceTrialIdx,:));
        
        flat_control_spatialFR = reshape(control_spatialFR.',1,[]);
        flat_dch_spatialFR = reshape(dch_spatialFR.',1,[]); 

        
        [autocor_cntrl,lags_cntrl] = xcorr(flat_control_spatialFR,200,'coeff');
        [autocor_dch,lags_dch] = xcorr(flat_dch_spatialFR,200,'coeff');
        
        MPP = 0.125;
%         figure(1); clf;
%         findpeaks(autocor_cntrl,'MinPeakProminence',MPP)
        [PKS_cntrl,~,~,prominence_cntrl] = findpeaks(autocor_cntrl,'MinPeakProminence',MPP);
        if isempty(prominence_cntrl)
            prominence_cntrl = 0;
        end
        
%         figure(2); clf;
%         findpeaks(autocor_dch,'MinPeakProminence',MPP)
        [PKS_dch,~,~,prominence_dch] = findpeaks(autocor_dch,'MinPeakProminence',MPP);
        if isempty(prominence_dch)
            prominence_dch = 0;
        end
        % find height of first autocorrelation peak 
        cells_autocorrelationScore(i) =  prominence_cntrl(1) - prominence_dch(1); 
    end
end