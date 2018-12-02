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

2. times_series.txt: 
each row of this file contains the time series of one voxel. The ordering of the voxels follows the same ordering of the file NoN_nodes_mod.txt. 

## Run file threshold_Cij.m
Input files: NoN_nodes_mod.txt and time_series.txt
Check that the $PATH for your input files is correct, if not, updated it. This code will create a series of correlation matrices thresholded at value lambda and save their plots.
## Run file NoN_construction_GC_two_lambdas
Input file: NoN_nodes_mod.txt
This code will output the final matrix J_NoN.txt which is the functional matrix of the network and its plot.

It will output a series of files:
- matrix_k_in.txt: list of in-degree for each module
- matrix_k_out.txt: list of out-degree for each pairs of modules
- matrix_outlinks.txt: list of number of out-links for each pairs of modules

## Run file circular_single_subj.m
Input file: NoN_nodes_mod.txt and matrix_outlinks.txt

This code will output a simple visualization of the network at the level of connections among modules. Each node is a module and a link between two nodes represents a weighted connection between two modules. 







