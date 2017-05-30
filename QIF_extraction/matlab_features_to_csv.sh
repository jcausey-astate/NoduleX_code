#!/bin/bash
#
# Converts output of lna.m results saved to a matlab file into 
# a csv file that can be further processed by other tools.
# The CSV will have one row per feature calculated, with the 
# first column as labels and additional columns containing values.

if [ $# -lt 2 ]
then
    echo "Usage: $0  matlab-file  csv-file  [field-delimiter]"
    echo "       matlab-file        name of matlab file to convert"
    echo "       csv-file           name of destination csv file"
    echo "       field-delimiter    character to separate fields (default=',')"
    exit 1
fi

delim="$3"
old_wd=$(pwd)
octave_wd=$(dirname "$0")

matfile_dir=$(dirname "$1")
csvfile_dir=$(dirname "$2")
matfile_base=$(basename "$1")
csvfile_base=$(basename "$2")

if [ "$matfile_dir" == "." ]
then
    matfile_dir="$old_wd"
fi

if [ "$csvfile_dir" == "." ]
then
    csvfile_dir="$old_wd"
fi

if [ -z "$3" ]
then
    delim=':'
fi

# Let octave read the matlab file and rip the cell array into individual CSV files
# per element (in case some elements contain more then one value)
cd "$octave_wd" && octave -Wi --eval "results_to_csv('${matfile_dir}/${matfile_base}', '${csvfile_dir}/${csvfile_base}', '${delim}')"
cd "$old_wd"
# That creates a large number of files with the pattern $2_part_NN.csv
# Collect them all into a single file named by $2
if [ -f "$2" ]
then
    rm "$2"
fi
for f in $2_part_*.csv 
do
    (cat "${f}"; echo) >> "$2" && rm "${f}"
done
# Collapse any multiple sets of values that were separated by multiple (2+) whitespaces
sed 's/   */,/g' "$2" > "$2.tmp" && mv "$2.tmp" "$2"
