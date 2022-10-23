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

#NWCHEM=/p/home/gkedz/Columbus/real_nwchem_for_columbus/bin/LINUX64/nwchem.dbg
export NWCHEM_TOP=/p/home/gkedz/Columbus/nwchem
export NWCHEM_TARGET=LINUX64
export NWCHEM_EXECUTABLE=/p/home/gkedz/Columbus/nwchem/bin/LINUX64/nwchem

mpiexec_mpt -n 2 $NWCHEM_EXECUTABLE khe_nwaoints.nw > khe_nwaoints.out
iwfmt.x > aoints.nw.sp
