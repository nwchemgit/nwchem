#!/usr/bin/env bash
if [[ $(uname -s) == "Darwin" ]]; then
    xcode_v=$(/usr/bin/xcodebuild -version |head -n1 |cut -d ' ' -f 2|cut -d '.' -f 1)
#    echo $xcode_v
    if [[ $( [ $xcode_v -ge 15 ] && echo 1) ]] ; then
	echo got xcode15
	export GOT_XCODE15=1
	export OMPI_FCFLAGS=" -Wl,-ld_classic "
        export OMPI_CFLAGS=" -Wl,-ld_classic -Wno-unused-command-line-argument "
#	export MPICH_FC="mpif90 -Wl,-ld_classic "
#	export MPICH_CC="mpicc -Wl,-ld_classic "
#	env|egrep MPICH_
    fi
fi
