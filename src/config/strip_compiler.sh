#/usr/bin/env bash
#      ifeq ($(shell basename -- $(FC)| cut -d \- -f 1),nvfortran)
#gcc  or gfortran?
if [[ "$1" == *fortran* ]] && [[ ! -z $_FC ]]; then
    echo $_FC
elif [[ "$1" == *cc* ]] && [[ ! -z $_CC ]]; then
    echo $_CC
elif [[ "$1" == *gfortran* ]] ; then
    echo gfortran
elif [[ "$1" == *flang* ]] ; then
    echo flang
else
    echo $(basename -- $1 | cut -d \- -f 1 | sed 's/[0-9]*//g')
fi

