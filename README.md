# Functional brain network construction from functional MRI signal
Construction of functional network from fMRI tasked-based data

You must have AFNI installed on your local machine.

## In AFNI:

Generate activation map by using 'correlation' as a way of thresholding:
If the name of the functional file was `functional.nii` (4D file, 3D coordinates
and 1D BOLD signal) this will generate a file  `functional.nii@+orig.BRIK` and
`functional.nii@+orig.HEAD`. These files
will be used in the following.

Extract the voxel coordinates from the activation map <br/>
Clusterize --> NN level = 2; Voxels = 30/40<br/>
-> Rpt -> Save the single modules one by one<br/>
Each module is save in a format
```
Clust_mask_0001+orig.BRIK
Clust_mask_0001+orig.HEAD
```
where 0001 indicates cluster 1.<br/>

### From terminal:

From terminal (awk): tranform the .BRIK/.HEAD file into a .txt file
```
1. $ 3dmaskdump Clust_mask_0001+orig. | awk '($4 !=0){print $0}' > cluster_1.txt
```
Output is in the format: x / y / z / module_label

Where module x, y, z are spatial coordinates of the voxel and module label is 1 for cluster 0001. This must be changed as needed, since the user might want to order the cluster at their convenience, i.e. the module label could be re-assigned as needed.

After saving all the clusters and converting them into .txt files, the user can proceed with cluster segmentation (if needed, not explained here).

All the segmented clusters are saved in a form:
```
clust_1.txt
clust_2.txt
.
.
.
clust_n.txt
```
Merge together all the clusters with
```
2. $ cat clust_1.txt clust_2.txt ... clust_n.txt > NoN_temp.txt
```
The NoN_temp.txt file contains the list of all the voxel coordinates and the corresponding module label, i.e x / y / z / module_label

Extract the value of each voxel's correlation value with the model used in the task.

```
3. $ 3dmaskdump functional.nii@+orig.[2] | awk -v th=VALUE '($4 > th || $4<-th){print $0}' > voxel_coord_th.txt
```

where VALUE is the value of the correlation threshold used in AFNI when generating the activation map

Join the voxel coordinates, with their correlation value and their module by using the bash script

```
4. $ ./nodes_coord_modules_and_corr.sh NoN_temp.txt voxel_coord_th.txt NoN_nodes_temp.txt
```

The output of the above command is NoN_nodes_temp.txt, structured as follows: x / y / z / correlation_val / module_label

Add the number of the line on the previous output, as first column:

```
5. $ awk '{print NR,$0}' NoN_nodes_temp.txt > NoN_nodes_mod.txt
```

NoN_nodes_mod.txt is a 6 columns file: NR / x / y / z / correlation value / module value, where NR = number of row.

Use the 4D file `functional.nii` and the bash script together with NoN_nodes_mod.txt to obtain the time series.

Get time series of the active voxels:

```
6. $ ./get_act_time_series.sh NoN_nodes_mod.txt letter_preprocessed.nii time_series.txt
```

The time_series.txt file contains the time series for each voxel in NoN_nodes_mod.txt

The files `NoN_nodes_mod.txt` and `time_series.txt` are further used to construct the functional network
associated the this functional MRI map. The structure of the files is

1. NoN_nodes_mod.txt :

6 columns file with NR / x / y / z / correlation value / module value


2. times_series.txt:
each row of this file contains the time series of one voxel. The ordering of the voxels follows the same ordering of the file NoN_nodes_mod.txt.

## Generate thresholded correlation matrices
Run file threshold_Cij.m

Input files: NoN_nodes_mod.txt and time_series.txt
Check that the $PATH for your input files is correct, if not, updated it. This code will create a series of correlation matrices thresholded at value lambda and save their plots.

## Generate Functional Network
Run file NoN_construction_GC_two_lambdas.m

Input file: NoN_nodes_mod.txt

This code will output the final matrix J_NoN.txt which is the functional matrix of the network and its plot.

It will output a series of files:
- matrix_k_in.txt: list of in-degree for each module
- matrix_k_out.txt: list of out-degree for each pairs of modules
- matrix_outlinks.txt: list of number of out-links for each pairs of modules

## Plot the functional network at the module level

Run file circular_single_subj.m

Input file: NoN_nodes_mod.txt and matrix_outlinks.txt

This code will output a simple visualization of the network at the level of connections among modules. Each node is a module and a link between two nodes represents a weighted connection between two modules.
