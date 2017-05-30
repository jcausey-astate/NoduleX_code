#!/bin/bash

function usage(){
	echo "Usage:"
	echo "    $(basename $0) [-g grey-dir] [-b bin-dir] [-o out-dir] nodule_list1 [nodule_list2 .. nodule_listN]"
	echo "        -g grey-dir   base directory for grey files (passed along to lna_batch_driver)"
	echo "        -b bin-dir    base directory for binary files (passed along to lna_batch_driver)"
	echo "        -o out-dir    base directory for output files (passed along to lna_batch_driver)"
	echo "        -l log-dir    base directory for log files"
	echo "        nodule_listX  text file containing a list of nodules to find features for; provide one"
	echo "                      such file for every parallel instance, so nodules from each file listed here"
	echo "                      will be processed in parallel with other lists."
	echo 
}

if [[ "$#" -lt 1 ]] ; then usage ; exit 1 ; fi

optionals=1

grey_dir=""
bin_dir=""
output_dir=""
log_dir="/dev/null"

while [ $optionals -eq 1 ] ; do
        if [ ${1:0:1} == "-" ] ; then
                param=${1:1:1}
                shift
                arg="$1"
                shift
                case $param in
                        g)
                          grey_dir="-g $arg"
                        ;;
                        b)
                          bin_dir="-b $arg"
                        ;;
                        o)
                          output_dir="-o $arg"
                        ;;
			l)
			  log_dir="$arg"
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

if [[ "$#" -lt 1 ]] ; then 
	echo "You must provide at least one nodule list file."
	usage 
	exit 2
fi

for i in $@ ; do
	run_log="${log_dir}/${i}_RUN.log"
	err_log="${log_dir}/${i}_ERR.log"
	if [ "$log_dir" == "/dev/null" ] ; then
		run_log="$log_dir"
		err_log="$log_dir"
	fi
	./lna_batch_driver.sh $grey_dir $bin_dir $output_dir `cat $i` >${run_log} 2>${err_log} &
done

echo "Started running $# lists in background at $(date )."

