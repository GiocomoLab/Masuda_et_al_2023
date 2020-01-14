function [I_sec, I_spike] = calculate_1DspatialInformation(singleCellSpatialFR,linearFractionalOccupancy)
    % Based on Skagg 1993
    % I_sec = ?(p_i * ?_i * log_2 (?_i/?))
    % I_spike = I_sec/?
    % p_i = linear Fractional Occupancy in spatial bin
    % ?_i = mean firing rate in the i-th spatial bin
    % ? = overall mean firing rate of the cell
    % 
    % Function Paramters
    % Inputs: 
    % - singleCellSpatialFR = trials x mean firing rate in each spatial bin
    % - linearFractionalOccupancy = occupancy in the ith bin/total recording time
    % Outputs:
    % - I_sec = bits/second
    % - I_spike = bits/spike

    if size(singleCellSpatialFR,2) ~= size(linearFractionalOccupancy,2)
       fprintf('Error: Number of spatial bins do not match between input');
       return
    end
    
    numBins = size(singleCellSpatialFR,2);
    
    p = linearFractionalOccupancy';
    lambda = nanmean(singleCellSpatialFR,1);
    lambda_total = nanmean(lambda);
    I_sec = 0;
    for i = 1:numBins
        I_sec = I_sec + (p(i) * lambda(i) * log2(lambda(i)/lambda_total));
        fprintf(num2str(I_sec));
    end
    
    I_spike = I_sec/lambda_total;
    
end