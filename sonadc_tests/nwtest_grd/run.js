#!/usr/bin/bash
#gen_pbs version: 1.2.0
#PBS -A ARLAP96070PET
#PBS -N nwqa
#PBS -j oe
#PBS -e nwqa.oe
#PBS -o nwqa.oe
#PBS -l walltime=0:10:00
#PBS -l select=1:ncpus=48:mpiprocs=48
#PBS -q debug

cd $PBS_O_WORKDIR

NWCHEM_EXECUTABLE=/p/home/gkedz/Columbus/nwchem/bin/LINUX64/nwchem
mpiexec_mpt -n 2 $NWCHEM_EXECUTABLE khe_nwgrd.nw > khe_nwgrd.out
