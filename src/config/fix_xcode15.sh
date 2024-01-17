#!/usr/bin/env bash
if [[ $(uname -s) == "Darwin" ]]; then
    xcode_v=$(clang --version|head -1 |cut -d ' ' -f 4 |cut -d . -f 1)
#    echo $xcode_v
    if [[ $( [ $xcode_v -ge 15 ] && echo 1) ]] ; then
	echo got xcode15
	export GOT_XCODE15=1
#figure out where the open-mpi libs are
        OMPILIBDIR=$(brew --cellar)/../opt/open-mpi/lib
	export OMPI_LDFLAGS=" -ld_classic  -L$OMPILIBDIR"
#	export OMPI_FCFLAGS=" -Wl,-ld_classic "
#        export OMPI_CFLAGS=" -Wl,-ld_classic -Wno-unused-command-line-argument "
#	export MPICH_FC="mpif90 -Wl,-ld_classic "
#	export MPICH_CC="mpicc -Wl,-ld_classic "
#	env|egrep MPICH_
    fi
fi
