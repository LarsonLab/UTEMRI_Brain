# flirt needs to be run in reference to another image 
# registration requires a reference

flirt -in /path/to/input/volume -ref /path/to/reference/volume -out /path/to/output/volume

 ## detail the nuances of reference matching