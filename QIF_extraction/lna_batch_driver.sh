#!/bin/bash

if [[ "$#" -lt 1 ]] ; then
	echo "Usage:"
	echo "    $(basename $0) -g grey-dir -b bin-dir -o out-dir  patient1 [patient2 .. patientN]"
	echo "        -g grey-dir            path to grey image files"
	echo "        -b bin-dir             path to binary image files"
	echo "        -o out-dir             path to output feature directory"
	echo "        patient1 .. patientN   Names of patients (or cases) that must correspond"
	echo "                               to subdirectories in the binary and grey directories"
	echo "                               under the main data directory (configure in script)."
	exit 1
fi

group_dir="analyze_format"
data_dir="$(pwd)"

data_dir="${data_dir}/${group_dir}"
bin_dir="${data_dir}/binary"
grey_dir="${data_dir}/grey"
output_dir="${data_dir}/features"

optionals=1

while [ $optionals -eq 1 ] ; do
	if [ ${1:0:1} == "-" ] ; then
		param=${1:1:1}
		shift
		arg="$1"
		shift
		case $param in
			g)
			  grey_dir="$arg"
			;;
			b)
			  bin_dir="$arg"
			;;
			o)
			  output_dir="$arg"
			;;
			*)
			  echo "Unknown option $param"
			  exit 1
			;;
		esac
	else
		optionals=0
	fi
done

if [[ ! -d "$output_dir" ]] ; then
	mkdir -p "$output_dir"
fi

for patient in $@ ; do
    nodule="$patient"
    patient="$(echo $patient | sed 's/\([^~]*\)~.*/\1/')"
    echo "Generating features for patient \"$patient\"..."
    patient_bin="${bin_dir}/${nodule}"
    patient_grey="${grey_dir}/${patient}"
    bin_file="$(ls $patient_bin/*.hdr)"
    grey_file="$(ls $patient_grey/*.hdr)"
    # Strip extensions:
    bin_file="${bin_file%.*}"
    grey_file="${grey_file%.*}"
    success=0
    octave -W -i --no-gui  --eval "lna_by_file('${bin_file}', '${grey_file}', '$nodule', '${output_dir}')" && (echo ; echo "OK" ; echo ; success=1 ) || (echo ; echo "FAILED" ; echo)
    if [ $success -eq 1 ]; then 
        rm "${patient_bin}/*.{diary,mat}"
    fi 
    echo "Finished patient \"$patient\"..."
done
