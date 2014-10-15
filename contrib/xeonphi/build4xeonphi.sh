#!/bin/bash -f
cd nwchem-6.5
export NWCHEM_TOP=`pwd`
export NWCHEM_TARGET=LINUX64
export USE_MPI=y
export USE_OPENMP=y
export USE_OFFLOAD=y
export NWCHEM_MODULES=tce
export BLASOPT=-mkl
export BLAS_SIZE=8
cd src
pwd
make nwchem_config  NWCHEM_MODULES="tce"
make FC=ifort 
