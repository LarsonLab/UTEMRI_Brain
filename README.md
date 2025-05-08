# Brain ultrashort-T2* component measurements
The purpose of this repoistory is to provide scripts and tools to generate fitted parameters to characterize the brain ultrashort-T2* component in healthy volunteers. It is currently implemented in MATLAB and includes code for reconstructing and fitting parameter maps from brain UTE scans.

It is an open-source collaborative platform to encourage anyone to contribute

## Citations

```
Boucneau T, Pao C, Tang S, Han M, Xu D, Henry RG, Larson PEZ. In vivo characterization of brain ultrashort-T2* components. Magnetic Resonance in Medicine. 2018 Aug;80(2):726-735. https://doi.org/10.1002/mrm.27037
```

```
Deveshwar N, Yao J, Han M, Dwork N, Shen X, Ljungberg E, Caverzasi E, Cao P, Henry R, Green A, Larson PEZ. Quantification of the in vivo brain ultrashort-T2* component in healthy volunteers. Magn Reson Med. 2024 Jun;91(6):2417-2430. doi: 10.1002/mrm.30013. Epub 2024 Jan 30. PMID: 38291598; PMCID: PMC11076235.
```

```
Shen X, Özen AC, Sunjar A, Ilbey S, Sawiak S, Shi R, Chiew M, Emir U. Ultra-short T2 components imaging of the whole brain using 3D dual-echo UTE MRI with rosette k-space pattern. Magn Reson Med. 2023 Feb;89(2):508-521. doi: 10.1002/mrm.29451. Epub 2022 Sep 25. PMID: 36161728; PMCID: PMC9712161.
```

```
Shen X, Caverzasi E, Yang Y, Liu X, Green A, Henry RG, Emir U, Larson PEZ. 3D balanced SSFP UTE MRI for multiple contrasts whole brain imaging. Magn Reson Med. 2024 Aug;92(2):702-714. doi: 10.1002/mrm.30093. Epub 2024 Mar 25. PMID: 38525680.
```


## Getting Started
1. Clone the `LarsonLab/UTEMRI_Brain` github repository
2. Download the sample dataset which includes the raw data (pfile)  and reconstructions of 18-degree, 12-degree, 6-degree scans respectively as well as a reconstructed UTE_AFI.mat file for the B1+ map correction


## Image Reconstruction
Run the script `recon_script` which requires access ot the raw Pfiles and a list of TEs for each flip angle (from zenodo). Image Reconstruction requires the Berkeley Advanced Reconstruction Toolbox (BART) to be installed. It can be found here: https://mrirecon.github.io/bart/

### Orchestra toolbox
This requires the Orchestra 1.7-1 sdk toolbox from GE which can be found here:https://www.gecares.com/s/GEHCForum. Extract the sdk into the same folder as the UTEMRI_Brain folder or add it to the matlab path with the command `addpath(genpath'/path/to/orchestra')`


## Generate ultrashort-T2* fitted parameter maps
Open `VFA_fitting_script.m` in MATLAB and follow through the code to run and generate fitted brain uT2 parameter maps. 

Approximate computation time on a 32-core CPU takes 1-1.5 hours
