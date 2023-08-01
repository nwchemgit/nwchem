#!/usr/bin/env bash
if [[ ! -z $GFORTRAN_MARCH ]] ; then
    exit=0
elif [[ -z $1 ]]; then
    #undefined
    exit=1
elif [[ $1 == "N" ]] || [[ $1 == "n" ]] || [[ $1 == "0" ]] || [[ ! -z $GFORTRAN_MARCH ]] ; then
    exit=0
else
    exit=1
fi
echo $exit
