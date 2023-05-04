# Masuda et al, 2023
Masuda, F.K., Sun, Y., Aery Jones, E.A., Giocomo, L.M. (2023, February). [Ketamine evoked disruption of entorhinal and hippocampal spatial maps](https://www.biorxiv.org/content/10.1101/2023.02.05.527227v1). bioRxiv doi:10.1101/2023.02.05.527227.

### System Requirements
* This code has been tested in MATLAB 2021b on Windows 10 and in MATLAB 2022b on Mac OS 13.3
* MATLAB codes requires the [gramm package](https://www.mathworks.com/matlabcentral/fileexchange/54465-gramm-complete-data-visualization-toolbox-ggplot2-r-like).

### Installation Guide
Estimated install time: 1 hour if MATLAB needs to be installed, 5 minutes otherwise.
1. Install [MATLAB](https://www.mathworks.com/help/install/install-products.html).
2. Clone this repository.
3. In MATLAB, add this repository and all its subfolders to your path:
`addpath(genpath('path\to\Masuda_et_al_2023'))`

### Instructions
1. Download preprocessed data from [Figshare](https://doi.org/10.6084/m9.figshare.22696309).
2. Run `masterScript.m` to combined data sessions into one struct. (Expected run time: 4 hours)
3. Generate plots from prepreprocessed MEC electrophysiology data using `plotAllCells.m` and scripts in _plottingFxns_ and _revisions_ directories. (Expected run time: a few minutes per script)
4. Generate plots from preprocessed hippocampal miniscope data using scripts in _miniscope_ (adapted from [the Sun et al, 2022 repository](https://github.com/yanjuns/Sun_Giocomo_2022_NComms))
