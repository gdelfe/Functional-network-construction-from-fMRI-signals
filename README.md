# NoN
Construction of functional network from fMRI tasked-based data

You must have AFNI installed on your local machine.

## From terminal:

extract the voxel coordinates from the activation map

1.$ 3dmaskdump activation_map_file[2] > voxel_coord.txt

Get rid of all the voxels below the correlation value used to generate the activatio map

2.$ awk -v th=VALUE '($4 > th || $4<-th){print $0}' voxel_coord.txt > voxel_coord_th.txt

Join the voxel coordinates, with their correlation value and their module

3.$ ./nodes_coord_modules_and_corr.sh NoN_temp.txt voxel_coord_th.txt NoN_nodes_temp.txt

Add number of line on the previous output

4.$ awk '{print NR,$0}' NoN_nodes_temp.txt > NoN_nodes_mod.txt

Convert BRICK and HEAD file to NIFTI

5.$ 3dAFNItoNIFTI -prefix letter_preprocessed.nii name_file_with_time_series_preprocessed

Get time series of the active voxels

6.$ ./get_act_time_series.sh NoN_nodes_mod.txt letter_preprocessed.nii time_series.txt

Output files further used to construct the functional network at the end of these above commands are:
1. NoN_nodes_mod.txt : 
6 columns file with NR / x / y / z / correlation value / module value
where NR - Number row
x, y, z - voxel coordinate
module value - value of the brain module for that voxel

2. times_series.txt: each row of this file contains the time series of one voxel. The ordering of the voxels follows the same ordering of the file NoN_nodes_mod.txt. 
