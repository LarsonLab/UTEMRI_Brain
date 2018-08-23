Moving from the MRI scanner to results

Using dcm2niix & FSL - BET, FAST, FLIRT 

BET+FAST for segmentation on T1 images, FLIRT (w/ BET) for registration between UTE and T1

An illustration of the pipeline that processes the images obtained from 3T MRI 

Working with volunteer data from T1 SPGR PRE-GAD (DCM folder 2)and UTE multiecho sat (DCM folder 3) to be compared against each other.

Results of dcm2nii conversion for both T1 and UTE are aggregated in preproc folder 

After BET/FAST processing (imaged through Papaya visualizer):
<img src="pics/post-fast.png">


Summer 2018: 
Registration for UTE_T2 data to T1 space. 
All data files can be found at /data/larson/brain_uT2/Rohit_Seg_Project. 
This also includes a sample first pass B/W segmentation using FSL. 
