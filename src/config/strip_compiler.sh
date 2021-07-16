#/usr/bin/env bash
#      ifeq ($(shell basename -- $(FC)| cut -d \- -f 1),nvfortran)
echo '###' > /tmp/dbg_`id -u`.txt
echo $(basename -- $1 | cut -d \- -f 1) >> /tmp/dbg_`id -u`.txt
echo $(basename -- $1 | cut -d \- -f 1)

