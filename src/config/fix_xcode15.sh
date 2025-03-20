#!/usr/bin/env bash
if [[ $(uname -s) == "Darwin" ]]; then
    xcode_v=$(clang --version 2>&1 |head -1 |cut -d ' ' -f 4 |cut -d . -f 1)
#    echo $xcode_v
    if [[ $( [ $xcode_v -ge 15 ] && echo 1) ]] ; then
	echo got xcode15
	export GOT_XCODE15=1
	#figure out where the open-mpi libs are
    if [ -x "$(command -v brew)" ]; then
        OMPILIBDIR=$(brew --cellar)/../opt/open-mpi/lib
	export OMPI_LDFLAGS=" -ld_classic  -L$OMPILIBDIR"
    fi
    fi
fi
