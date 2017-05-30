#!/bin/bash
# Converts results from lna feature production to individual CSV files.

usage="Usage:\n\t  $(basename ${0}) [--output_dir output_dir] matlab-file [matlab-file ...]"

if [ $# -lt 2 ]
then
    echo -e $usage
    exit 1
fi

output_dir=""

if [ "$1" == "--output_dir" ]
then
    output_dir="${2}/"
    shift
    shift
fi

source_dir=$(dirname $0)

for file in $@ 
do
    outname=$(basename $file)
    outname=${outname%.*}
    outname="${outname}.csv"
    bash ${source_dir}/matlab_features_to_csv.sh "$file" "${output_dir}$outname" >/dev/null || exit 2
done
