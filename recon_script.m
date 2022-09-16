%  UTE Brain Reconstruction script
% Input: Pfile and TE list
% Output: Two 4-D reconstructions (height, width, num_slices, num_echoes)

% try 440 Hz offset?

addpath(genpath('/home/sf719662/Documents/UTEMRI_Brain'));

% Define pfile list and corresponding set of TEs to use in reconstruction
pfilename = {'P38400.7_11211332', 'P39424.7_11211346', 'P40960.7_11211359'};
TE_set = {'multi_utes-nd_8TE_18-deg.dat', 'multi_utes-nd_8TE_12-deg.dat', 'multi_utes-nd_8TE_06-deg.dat'};
%% reconstruct images
% 300 Hz offset
offset_frequency = 300;

for n = 1:length(pfilename)
    ute_brain_recon(pfilename{n}, TE_set{n}, offset_frequency);
end


