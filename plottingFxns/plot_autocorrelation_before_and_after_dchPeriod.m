function plot_autocorrelation_before_and_after_dchPeriod(cells)

plotWidth = 500;
plotHeight = 1000;
figure('Position',[200 200 plotWidth plotHeight]); hold on;
for i = 1:numel(cells.metadata)
    clf;
    smoothFR = cells.spatialFRsmooth;
    
    dch = cells.dch(i).dch;
    decoherenceTrialIdx = dch.decoherenceIdx;
    
    if ~isempty(decoherenceTrialIdx)
        dchIndxLength = numel(decoherenceTrialIdx);
%         controlTrials = 101-dchIndxLength:100;
        controlTrials = 1:dchIndxLength;
        if controlTrials(1) < 1
            controlTrials = 1:dchIndxLength;
        end
        control_spatialFR = squeeze(smoothFR(i,controlTrials,:));
        dch_spatialFR = squeeze(smoothFR(i,decoherenceTrialIdx,:));
        
%         norm_control_spatialFR = normalize(control_spatialFR,2);
%         norm_dch_spatialFR = normalize(dch_spatialFR,2);
%         
        flat_control_spatialFR = reshape(control_spatialFR.',1,[]);
        flat_dch_spatialFR = reshape(dch_spatialFR.',1,[]); 

        flat_control_spatialFR = fillmissing(flat_control_spatialFR,'previous');
        flat_dch_spatialFR = fillmissing(flat_dch_spatialFR,'previous');
        
        [autocor_cntrl,lags_cntrl] = xcorr(flat_control_spatialFR,400,'coeff');
        [autocor_dch,lags_dch] = xcorr(flat_dch_spatialFR,400,'coeff');
       
        %%
        MPP = 0.125;
        [~,~,~,prominence_cntrl] = findpeaks(autocor_cntrl,'MinPeakProminence',MPP,'SortStr','descend');
        if isempty(prominence_cntrl)
            prominence_cntrl = 0;
        end
        
        %%
        [~,~,~,prominence_dch] = findpeaks(autocor_dch,'MinPeakProminence',MPP,'SortStr','descend');
        if isempty(prominence_dch)
            prominence_dch = 0;
        end
        %
        if numel(prominence_dch)<2 && numel(prominence_cntrl)<2
            cells_autocorrelationScore = 0;
            p_cntrl = 0;
            p_dch = 0;
        elseif numel(prominence_dch)<2
            % find height of second autocorrelation peak 
            cells_autocorrelationScore =  prominence_cntrl(2); 
            p_dch = 0;
        elseif numel(prominence_cntrl)<2
            cells_autocorrelationScore =  0 - prominence_dch(2); 
            p_cntrl = 0;
        else
             % find height of second autocorrelation peak 
            
            p_cntrl = prominence_cntrl(2);
            p_dch = prominence_dch(2);
            cells_autocorrelationScore =  p_cntrl-p_dch/p_cntrl; 
        end
        %%
        clf;
        t = tiledlayout(3,1);
        colormap('hot')
        title(t,sprintf('Control vs. Decoherence Trials Cell %i',i))
        nexttile();
        imagesc(control_spatialFR)
        title('Control');
 
        nexttile();
        imagesc(dch_spatialFR)
        title('Decoherence');
         
        nexttile();
        plot(lags_cntrl(401:end),autocor_cntrl(401:end));
        hold on;
        plot(lags_dch(401:end),autocor_dch(401:end));
        title(sprintf('Control(%0.2f) vs Decoherence(%0.2f) Autocorrelation: %0.2f',p_cntrl,p_dch,cells_autocorrelationScore))
        legend('Control','Decoherence')
        
    else
        fprintf('No Decoherence period\n') 
    end
    pause
end