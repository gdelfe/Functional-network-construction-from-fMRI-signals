# bash script -- Gino Del Ferraro -- January 2017
#
# Argument 1: file with the coordinate of the active voxels in format:
#  NR/ i / j / k / correlation value / module
#
# Argument 2: file NIFTI 4D containing the functional image
#
# Argument 3: output file

if [ $# -eq 3 ]; then # exactly 3 argument was passed..use it..its available in $1, $2 and $3
    echo "Arguments: $1 $2 $3"
else # either 0 or <3 arguments were passed...error out.
echo "Incorrect number of arguments passed"
exit 1
fi


rm "$3"
while read -r line # read the file in argument 1, i.e activated voxel coordinate
do
    i=$(echo $line | awk '{print $2}') # read i
    j=$(echo $line | awk '{print $3}') # read j
    k=$(echo $line | awk '{print $4}') # read k
nifti_tool -disp_ts $i $j $k -input "$2" | tail -n 1 >> "$3" # print the corresponding time series of i,j,k into a file
done < "$1"


