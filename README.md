# NoN
Construction of functional network from fMRI tasked-based data

You must have AFNI installed on your local machine.

From terminal:
# extract the voxel coordinates from the activation map
1. 3dmaskdump activation_map_file[2] > voxel_coord.txt
2. awk -v th=VALUE '($4 > th || $4<-th){print $0}' voxel_coord.txt > voxel_coord_th.txt

