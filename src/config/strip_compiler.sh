#/usr/bin/env bash
#      ifeq ($(shell basename -- $(FC)| cut -d \- -f 1),nvfortran)
echo $(basename -- $1 | cut -d \- -f 1 | sed 's/[0-9]*//g')

