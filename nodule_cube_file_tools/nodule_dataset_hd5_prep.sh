#!/bin/bash

# Prepare HD5 nodule file by assigning classes according to 
# malignancy score and sorting so that the (+) class comes first.

set -e

nodule_file_tools=$dirname($0)

if [[ $# -lt 7 ]] ; then
	echo "usage: $(basename $0) nodule-file.hd5 m0 m1 m2 m3 m4 m5"
	echo ""
	echo "    You must provide hd5 file and 6 malignancy mapping scores, in order"
	echo "    as follows: m0 m1 m2 m3 m4 m5 , where each score is the class to"
	echo "    assign to any nodule matching that malignancy score.  Use -1 to"
	echo "    mark a malignancy level as excluded."
	exit 1
fi

FILE="$1"

python ${nodule_file_tools}/get_set_nodule_classes.py -m $2 $3 $4 $5 $6 $7  ${FILE}
python ${nodule_file_tools}/sort_nodules_in_hd5.py -a ids ${FILE} 
python ${nodule_file_tools}/sort_nodules_in_hd5.py -a classes -r ${FILE}  

echo "Nodule file ${FILE} prepared OK"

