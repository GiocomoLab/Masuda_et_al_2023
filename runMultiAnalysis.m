% sessions = {...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190619_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190620_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190621_keicontrasttrack_ketamine2_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190623_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/G4/G4_190624_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCN1/HCN1_190618_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCN1/HCN1_190619_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCN1/HCN1_190620_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCN1/HCN1_190621_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCN1/HCN1_190623_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/B1/B1_190527_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/B1/B1_190529_keicontrasttrack_ketamine1_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/B1/B1_190528_keicontrasttrack_ketamine2_g0',...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/B3/B3_0515_contrasttrack_ketamine1', ...
% '/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/B3/B3_0516_contrasttrack_ketamine1', ...
%         };

sessions = {'/Volumes/groups/giocomo/export/data/Projects/JohnKei_NPH3/HCN1/HCN1_190625_keicontrasttrack_propofol1_g0'};
for n = 1:numel(sessions)
    try
        singleSessionRasterplots(sessions{n},session_name, trackLength, plotWidth, plotHeight, preDrugTrials)
        drawLicksMultiSessions(sessions{n})
        
    catch
        warning(sprintf('FAILED: %s\n',sessions{n}))
    end
end

