#!/usr/bin/env bash
# check if environment USE_PEIGS is
# compatible with GA generated _USE_PEIGS
NWCHEM_TOP=$1
file_check=$NWCHEM_TOP/src/peigs_check_done.txt
if [ -f $file_check ]; then
    cat $file_check
    exit 0
fi
if [[  -z "${USE_PEIGS}" ]]; then
    echo 0 >  $file_check
    exit 0
fi
#extract _USE_PEIGS from ga_use_peigs.txt
if [ -f $NWCHEM_TOP/src/ga_use_peigs.txt ]; then
    _USE_PEIGS=$(cat $NWCHEM_TOP/src/ga_use_peigs.txt)
else
    exit 0
fi
_use_peigs_set=n
if [[ "$USE_PEIGS" =~ ^(y|Y|1)$ ]]; then
    use_peigs_set=Y
else
    use_peigs_set=N
fi
#echo use_peigs_set "${use_peigs_set}"
if [[ "${use_peigs_set}" != "${_USE_PEIGS}" ]]; then
    echo 1 >  $file_check
fi
# 
