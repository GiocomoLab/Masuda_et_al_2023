function interneuronFlag = findInterneurons(cells)
% Add interneuron flag (=
% 1 if baseline trials 1:50 has a mean FR > 15hz
nCells = size(cells.spatialFRsmooth,1);
interneuronFlag = nan(nCells,1);
baselineFRsmooth = cells.spatialFRsmooth(:,1:50,:);
baselineFRsmooth = permute(baselineFRsmooth, [3 2 1]);
linearized_cellFR = squeeze(reshape(baselineFRsmooth,[],1,nCells));

mean_linearized_cellFR = mean(linearized_cellFR,1);
interneuronFlag = mean_linearized_cellFR>15;
interneuronFlag = interneuronFlag';
                                              
