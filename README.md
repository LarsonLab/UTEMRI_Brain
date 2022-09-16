# Brain ultrashort-T2 component measurements
The purpose of this repoistory is to provide scripts and tools to generate fitted parameters to characterize the brain ultrashort-T2 component in healthy volunteers. It is currently implemented in MATLAB and includes code for reconstructing(?) and fitting parameter maps from brain UTE scans.

It is an open-source collaborative platform to encourage anyone to contribute

## Citations
Boucneau T, Pao C, Tang S, Han M, Xu D, Henry RG, Larson PEZ. In vivo characterization of brain ultrashort-T2 components. Magnetic Resonance in Medicine. 2018 Aug;80(2):726-735. https://doi.org/10.1002/mrm.27037


## Getting Started
1. Clone the `LarsonLab/UTEMRI_Brain` github repository
2. Download the sample dataset which includes the raw data (pfile)  and reconstructions of 18-degree, 12-degree, 6-degree scans respectively as well as a reconstructed UTE_AFI.mat file for the B1+ map correction

### Orchestra toolbox
This requires the Orchestra 1.7-1 sdk toolbox from GE which can be found here:https://www.gecares.com/s/GEHCForum. Extract the sdk into the same folder as the UTEMRI_Brain folder or add it to the matlab path with the command `addpath(genpath'/path/to/orchestra')`


## Example fitting script
Open `VFA_fitting.m' script in MATLAB and follow through the code to run and generate fitted brain uT2 parameter maps. 

Approximate computation time on a 32-core CPU takes 1-1.5 hours
