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

NWCHEM=/p/home/gkedz/Columbus/real_nwchem_for_columbus/bin/LINUX64/nwchem.dbg
mpirun -np 1 $NWCHEM KHe_70_nwaoints.nw > KHe_70_nwaoints.out
iwfmt.x > aoints.nw.sp
