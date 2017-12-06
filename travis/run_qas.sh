#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 nprocs=2
 os=`uname`
 if [[ "$ARMCI_NETWORK" == "MPI-PR" ]]; then
     nprocs=3
    if [[ "$os" == "Darwin" && "$MPI_IMPL" == "openmpi" ]]; then 
	export MPIRUN_NPOPT="-mca mpi_yield_when_idle 0 --oversubscribe -np "
    fi
 fi
 if [[ "$os" == "Linux" && "$MPI_IMPL" == "mpich" ]]; then
    export MPIRUN_PATH=/usr/bin/mpirun.mpich
 fi
 if [[ "$NWCHEM_MODULES" == "tce" ]]; then
   cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs tce_n2 tce_ccsd_t_h2o tce_h2o_eomcc
   cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs tce_ipccsd_f2 tce_eaccsd_ozone
 else
     cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs $nprocs dft_he2+ prop_mep_gcube
     cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs cosmo_h2o_dft  h2o_opt
     cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs dft_siosi3
     cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs pspw 
     if [[ "$ARMCI_NETWORK" != "MPI-PT" ]] || [[ "$os" == "Darwin" && "$ARMCI_NETWORK" != "MPI-PR" ]]; then
       cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs tddft_h2o h2o2-response
       cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs pspw_md
     fi
 fi
