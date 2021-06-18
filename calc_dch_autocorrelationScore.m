function [cells_autocorrelationScore,p_cntrl, p_dch] = calc_dch_autocorrelationScore(cells)

cells_autocorrelationScore = nan(size(cells.metadata,1),1);
p_cntrl = nan(size(cells.metadata,1),1);
p_dch = nan(size(cells.metadata,1),1);
for i = 1:size(cells.metadata,1)
   
    smoothFR = cells.spatialFRsmooth;
    try
        dch = cells.dch(i).dch;
        decoherenceTrialIdx = dch.decoherenceIdx;
    end
    
    if ~isempty(decoherenceTrialIdx)
        dchIndxLength = numel(decoherenceTrialIdx);
        controlTrials = 101-dchIndxLength:100;
%         controlTrials = 1:dchIndxLength;
        if controlTrials(1) < 1
            controlTrials = 1:dchIndxLength;
        end
        control_spatialFR = squeeze(smoothFR(i,controlTrials,:));
        dch_spatialFR = squeeze(smoothFR(i,decoherenceTrialIdx,:));
        
        
        
        norm_control_spatialFR = normalize(control_spatialFR,2);
        norm_dch_spatialFR = normalize(dch_spatialFR,2);
        
        flat_control_spatialFR = reshape(norm_control_spatialFR.',1,[]);
        flat_dch_spatialFR = reshape(norm_dch_spatialFR.',1,[]); 

        flat_control_spatialFR = fillmissing(flat_control_spatialFR,'previous');
        flat_dch_spatialFR = fillmissing(flat_dch_spatialFR,'previous');
        
        [autocor_cntrl,lags_cntrl] = xcorr(flat_control_spatialFR,200,'coeff');
        [autocor_dch,lags_dch] = xcorr(flat_dch_spatialFR,200,'coeff');
        
%         figure(1); clf; hold on;
%         plot(lags_cntrl,autocor_cntrl)
%         plot(lags_dch,autocor_dch)
        
        MPP = 0.125;
%         figure(1); clf;
%         findpeaks(autocor_cntrl,'MinPeakProminence',MPP)
        [PKS_cntrl,~,~,prominence_cntrl] = findpeaks(autocor_cntrl,lags_cntrl,'MinPeakProminence',MPP,'SortStr','descend');
        if isempty(prominence_cntrl)
            prominence_cntrl = 0;
        end

%         figure(2); clf;
%         findpeaks(autocor_dch,'MinPeakProminence',MPP)
        [PKS_dch,~,~,prominence_dch] = findpeaks(autocor_dch,lags_dch,'MinPeakProminence',MPP,'SortStr','descend');
        if isempty(prominence_dch)
            prominence_dch = 0;
        end
        %%
        if numel(prominence_dch)<2 && numel(prominence_cntrl)<2
            cells_autocorrelationScore(i) = nan;
            p_cntrl(i) = nan;
            p_dch(i) = nan;
        elseif numel(prominence_dch)<2
            % find height of second autocorrelation peak 
            cells_autocorrelationScore(i) =  prominence_cntrl(2); 
            p_dch(i) = nan;
        elseif numel(prominence_cntrl)<2
            cells_autocorrelationScore(i) =  0 - prominence_dch(2); 
            p_cntrl(i) = nan;
        else
             % find height of second autocorrelation peak 
            
            p_cntrl(i) = prominence_cntrl(2);
            p_dch(i) = prominence_dch(2);
            cells_autocorrelationScore(i) =  (p_cntrl(i)-p_dch(i))/(p_cntrl(i)); 
        end
    end
end