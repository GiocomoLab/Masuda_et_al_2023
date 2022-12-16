addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/MalcolmFxn/'));
addpath(genpath('/Volumes/groups/giocomo/export/data/Users/KMasuda/Neuropixels/npy-matlab'));

% myKsDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCN1/HCN1_190621_keicontrasttrack_ketamine1_g0/HCN1_190621_keicontrasttrack_ketamine1_g0_imec0/';
% myKsDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190619_keicontrasttrack_ketamine1_g0/G4_190620_keicontrasttrack_ketamine1_g0_imec0/';
% myKsDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G5/G5_190705_keicontrasttrack_ketamine1_g0/G5_190705_keicontrasttrack_ketamine1_g0_imec0';
myKsDir = '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCNe1/HCNe1_190830_keicontrasttrack_ketamine1_g0/HCNe1_190830_keicontrasttrack_ketamine1_g0_imec0';
lfpD = dir(fullfile(myKsDir, '*.lf.bin')); % LFP file from spikeGLX specifically
lfpFilename = fullfile(myKsDir, lfpD(1).name);


lfpFs = 2500;  % neuropixels phase3a
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX

[lfpByChannel, allPowerEst, F, allPowerVar] = ...
    lfpBandPowerTime(lfpFilename, lfpFs, nChansInFile, []);

% [lfpByChannel, allPowerEst, F, allPowerVar] = ...
%     lfpBandPower(lfpFilename, lfpFs, nChansInFile, []);

chanMap = readNPY(fullfile(myKsDir, 'channel_map.npy'));
nC = length(chanMap);

allPowerEst = allPowerEst(:,chanMap+1)'; % now nChans x nFreq

% plot LFP power
dispRange = [0 80]; % Hz
marginalChans = [10:50:nC];
freqBands = {[1.5 4], [4 10], [10 30], [30 80], [80 200]};

plotLFPpower(F, allPowerEst, dispRange, marginalChans, freqBands);