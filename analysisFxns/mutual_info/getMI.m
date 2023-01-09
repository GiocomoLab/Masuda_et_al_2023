function [MI_XY, p_X, p_Y] = getMI(X, Y, X_edges, Y_edges)
% Determine the mutual information between spikes and position.
% Uses only spikes for speed > 2cm/s
%
% Input variables should be in the form of: 
%       X in counts over time, eg. #spikes/sec
%       Y in value over time, eg. position at each time-point
%       *NOTE: X and Y should be aligned and the same length*
%
%       X_edges: vector defining min:bin:max for X histogram
%       Y_edges: vector defining min:bin:max for Y histogram
%       num_bins: number of time bins used to bin X,Y
%
% Output entropy_X = the general variability in X
%        entropyX_givenY = the amount of variability unexplained by Y
%        MI_XY = mutual information between X and Y
%
% Created by Isabel Low 4/16/18
% Last edited by IL 4/20/18

close all
num_bins = numel(X);

% p(x,y): joint probability of x and y
p_XY = histcounts2(X,Y,X_edges,Y_edges)/num_bins;
p_X = sum(p_XY,2);
p_Y = sum(p_XY,1);

% determine MI
p_X_times_p_Y = repmat(p_X, 1, length(p_Y)) .* repmat(p_Y, length(p_X), 1);
MI_mat = p_XY .* log2(p_XY ./ p_X_times_p_Y);
MI_XY = nansum(MI_mat(:));