# JohnKeiNPAnalysis
Neuropixel Analysis code specific to Kei and John

## File Structure
calcSpeed.m = universalFXN
drawSpeed.m = universalFXN **Change OAK path**

calcFRmapCorrMatrixAllCells.m 
- calcSmoothedFR_SpatialBin.m


runMultiAnalysis.m **comment this**

runMultiPostProcessing.m **comment this + FXN + FLAG FOR STITCH + SAVE MI VECTORS**
	- stitchSynchedNPdata.m
		- concatenateNPMatFiles.m
			- concatenateSpStructs.m
	- singleSessionRasterplots.m ** add save metadata (depths, waveform, cells to plot)
	- combineRASTERS.m
	- drawLicksSingleSessions.m = universalFXN **Change OAK path**
	- getMIAllCells.m = unfold into parent fxn
		- getMI.m
		- calcSmoothedFR_Time.m 



gauss_smoothing.m =  **delete from master**

plotFiringRate.m **need to transform into a plotAvgFR for all sessions using calcSmoothedFR_Time.m**

singleSessionRasterplotsNaNspikes.m = DELETE 
drawLicksMultiSessions.m = DELETE

baselineCntrlKetamine_rasterPlots.m = indivdual script 
dosageRasters.m = individual script

scratch.m = individual script
scratch_for_fields.m = individual script
scratch_for_repeating.m = individual script
singleSessionRasterplotsSelectTrials.m = individual script