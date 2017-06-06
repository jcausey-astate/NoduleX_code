#!/usr/bin/env bash
#
# Walk LIDC-IDRI directory structure looking for directories with more than N dicom files in it
# Echo any such directory path to stdout.
# N will default to 20, but can be supplied as second arg.
# If you set N to a negative value, you will find directories with LESS THAN |N| dicom files, but
# containing AT LEAST ONE regular file.
# 
# usage:
# walk-doi.sh  root_dir  N

if [ "$1" = "-h" ]
then
    echo "Usage:  $(basename $0) ROOT-DIR N"
    echo 
    echo "        ROOT-DIR:  Sample's top-level directory."
    echo "        N:         Minimum number of dicom files required to include directory"
    echo "                   or maximum number if set to negative (ex -20 finds"
    echo "                   directories with less than 20 dicom files, and at least ONE"
    echo "                   regular file of any kind."
    echo
    exit 1
fi

# Platform-independant substitute for "readlink -f" based on code found at:
# http://stackoverflow.com/a/19250873
# Modified to only care about "existence" (file or directory is OK)
function get_realpath() {

    [[ ! -e "$1" ]] && return 1 # failure : file does not exist.
    [[ -n "$no_symlinks" ]] && local pwdp='pwd -P' || local pwdp='pwd' # do symlinks.
    echo "$( cd "$( echo "${1%/*}" )" 2>/dev/null; $pwdp )"/"${1##*/}" # echo result.
    return 0 # success

}

ROOTDIR=$1
N=20

if [ "$#" -gt 1 ] ; then
    N=$2
fi

if [ "$#" -lt 1 ] ; then
    ROOTDIR=$(pwd)
fi

LESS_THAN=0

if [ "$N" -lt "0" ] ; then
    LESS_THAN=1
    N=$(($N * -1))
fi

#DIR=$(readlink -f $ROOTDIR)  # NOTE: Will not work on Mac
DIR=$(get_realpath $ROOTDIR)  #       This one is platform-agnostic
LEN=${#DIR}
LEN=$(( LEN - 1 ))
if [[ "${DIR:$LEN}" == "/" ]] ; then
    DIR=${DIR:0:$LEN}             # drop the extra trailing slash
fi

for dir in `find $DIR -type d` ; do
    NDICOMS=$(ls -l $dir/*.{dcm,dicom} 2> /dev/null | wc -l)
    if [ "$LESS_THAN" -eq 0 ] ; then
        if [ "$NDICOMS" -gt "$N" ] ; then
            echo "$dir"
        fi
    else
        NDIRS=$(ls -ld $dir/*/ 2> /dev/null | wc -l)
        NFILES=$(ls -l $dir 2> /dev/null | wc -l)
        NFILES=$(($NFILES - 1)) # Account for "total" line
        if [ "$NDICOMS" -lt "$N" ] && [ "$NFILES" -gt 0 ] && [ "$NDIRS" -ne "$NFILES" ] && [ "$dir" != "$DIR" ] ; then
            echo "$dir"
        fi
    fi
done
