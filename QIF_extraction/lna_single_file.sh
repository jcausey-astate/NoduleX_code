#!/bin/bash

# Usage:
#    lna_single_file.sh grey_dir binary_dir patient_id output_dir

if [ $# -lt 1 ] ; then
	echo "Usage:"
	echo "     $(basename ${0}) grey_dir binary_dir patient_id output_dir"
	echo
	exit 1 
fi

bin_dir="$2"
grey_dir="$1"
output_dir="$4"
patient="$3"

echo "Generating features for patient \"$patient\"..."
    
octave -W -i --no-gui  --eval "lna_by_file('${bin_dir}', '${grey_dir}', '$patient', '${output_dir}')" && (echo ; echo "OK" ; echo ) || (echo ; echo "FAILED" ; echo)
    
echo "Finished patient \"$patient\"..."

