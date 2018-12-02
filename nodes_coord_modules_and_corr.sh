# This script merges the nodes which are in the clusters with their coordinates.
#
# bash script -- Gino Del Ferraro -- May 2017
#
# Argument 1: file with the mask: x / y / z / module
#
# Argument 2: file with the voxel coordinates (scan) of the act voxel above the threshold: x / y / z / correlation
#
# Argument 3: output file

if [ $# -eq 3 ]; then # exactly 3 argument was passed..use it..its available in $1, $2 and $3
echo "Arguments: $1 $2 $3"
else # either 0 or <3 arguments were passed...error out.
echo "Incorrect number of arguments passed"
exit 1
fi


rm "$3"
while read -r line # read the file in argument 1, i.e node in the module
do
i=$(echo $line | awk '{print $1}') # read i
j=$(echo $line | awk '{print $2}') # read j
k=$(echo $line | awk '{print $3}') # read k
module=$(echo $line | awk '{print $4}') # read module
grep  "$i $j $k" "$2" | awk '{print $0,'$module'}' >> "$3" # print coordinate x y z, correlation, module
done < "$1"
