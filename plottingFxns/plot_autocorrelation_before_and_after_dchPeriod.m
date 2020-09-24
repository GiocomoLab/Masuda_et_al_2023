function plot_autocorrelation_before_and_after_dchPeriod(cells)


addpath(genpath('/Users/KeiMasuda/Documents/MATLAB/Add-Ons/Functions/gramm (complete data visualization toolbox, ggplot2_R-like)/code'));
%%

for i = 256:numel(cells.metadata)
    clf;
    smoothFR = cells.spatialFRsmooth;
    
    dch = cells.dch(i).dch;
    decoherenceTrialIdx = dch.decoherenceIdx;
    
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
       
%         controlStability = calculateStabilityScore(control_spatialFR);
%         dchStability = calculateStabilityScore(dch_spatialFR);
        
        t = tiledlayout(2,2);
        colormap('hot')
        title(t,sprintf('Control vs. Decoherence Trials Cell %i',i))
        nexttile();
        imagesc(control_spatialFR)
        title('Control');
 
        nexttile();
        imagesc(dch_spatialFR)
        title('Decoherence');
         
        nexttile([1 2]);
        plot(lags_cntrl,autocor_cntrl);
        hold on;
        plot(lags_dch,autocor_dch);
        title('Control vs Decoherence Autocorrelation')
        
%         nexttile();
%         plot(lags_dch,autocor_dch);
%         title('Dechoerence Autocorrelation')
    else
        fprintf('No Decoherence period\n') 
    end
    pause
end