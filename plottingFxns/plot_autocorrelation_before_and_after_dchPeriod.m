function plot_autocorrelation_before_and_after_dchPeriod(cells)

% function constants 
plotWidth = 500;
plotHeight = 1000;
spatialBinNum = 200;
MPP = 0.125; %Minimum Peak Prominence

numCells = size(cells.metadata,1);
drawFig = true;


% declare global vaariables
all_pcntrl = NaN(numCells,1);
all_pdch = NaN(numCells,1);
all_cells_autocorrelationScore = NaN(numCells,1);

% if drawFig
%     figure('Position',[200 200 plotWidth plotHeight]); hold on;
% end
%%
for i = 1:numCells

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
        
        %% Calculate autocorrelations 
         [autocor_cntrl,lags_cntrl] = xcorr(flat_control_spatialFR,400,'coeff');
        [autocor_dch,lags_dch] = xcorr(flat_dch_spatialFR,400,'coeff');
        %%
        % Find Peaks in autocorrelation curves
        [pks_cntrl,~,lags_cntrl,prominence_cntrl] = findpeaks(autocor_cntrl,'MinPeakProminence',MPP,'SortStr','descend');
        autocor_cntrl = autocor_cntrl(1:numel(autocor_cntrl));
        if isempty(prominence_cntrl)
            prominence_cntrl = 0;
        end
        
        %
        [pks_dch,~,lags_dch,prominence_dch] = findpeaks(autocor_dch,'MinPeakProminence',MPP,'SortStr','descend');
        autocor_dch = autocor_dch(1:numel(autocor_dch));
        if isempty(prominence_dch)
            prominence_dch = 0;
        end
        %% calculate 2nd peak (not the self correlation peak)
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
             % find dF/F for second autocorrelation peak peak 
            p_cntrl = prominence_cntrl(2);
            p_dch = prominence_dch(2);
            cells_autocorrelationScore =  p_cntrl-p_dch/p_cntrl; 
        end
        
        %%%%%%%%%%%%%
        %% LBQ TEST 
        %%%%%%%%%%%%%
%         residuals_cntrl = flat_control_spatialFR - mean(flat_control_spatialFR);
%         [h,pValue,stat,cValue] = lbqtest(residuals_cntrl,'lags',200);
% %         archtest(residuals_cntrl)
% %         morans_I(control_spatialFR,ones(size(control_spatialFR)))
% 
%         residuals_dch = flat_dch_spatialFR - mean(flat_dch_spatialFR);
%         [h,pValue,stat,cValue]=lbqtest(residuals_dch,'lags',200);
% %         archtest(residuals_dch)
        %%
        if drawFig
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
            plot(autocor_cntrl);
            hold on;
            plot(autocor_dch);
            title(sprintf('Control(%0.2f) vs Decoherence(%0.2f) Autocorrelation: %0.2f',p_cntrl,p_dch,cells_autocorrelationScore))
            legend('Control','Decoherence')
            pause;
        end
        %%
        all_pcntrl(i) = p_cntrl;
        all_pdch(i) = p_dch;
        all_cells_autocorrelationScore(i) = cells_autocorrelationScore;
        fprintf('finished %i iterations\n', i);   
    else
        fprintf('No Decoherence period\n') 
    end
end


%% Summary Plots

[h,p] = ttest2(all_pcntrl, all_pdch);
fprintf('The hypothesis test is %i; P-value is %.02d\n',h,p);
% Compare Before Dch and After Dch Figures
close all; 
data = vertcat(all_pcntrl,all_pdch);
nanIndx = isnan(data);
data = data(~nanIndx);

metadata = cells.metadata(:,2);
metadata = vertcat(metadata,metadata);
metadata = metadata(~nanIndx);

color = vertcat(...
    repmat({'Control'}, numCells),...
    repmat({'Decoherence'}, numCells)...
    );
color = color(~nanIndx);
%
clear g;
g=gramm('x',color,'y',data,'color',color);
% g.stat_bin();
% g.stat_violin('normalization','width','dodge',0,'fill','edge','half','true');
g.stat_boxplot('notch','true','geom_jitter');
% g.stat_summary();
g.set_title('Autocorrelation Score Control vs Decoherence Period');
g.set_names('x','','y', 'Hz');
customColorMap = [ % 0.5 0.5 0.5 %grey
    0.8 0.2 0.8 %magenta
    0 0.8 0.2]; %green
g.set_color_options('map',customColorMap);
g.draw()