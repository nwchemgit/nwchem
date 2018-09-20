#!/bin/bash -f
# Exit on error
set -ev
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 os=`uname`
 nprocs=2
 do_largeqas=1
 if [[ "$os" == "Linux" && "$MPI_IMPL" == "mpich" ]]; then
    export MPIRUN_PATH=/usr/bin/mpirun.mpich
 fi
 if [[ "$os" == "Darwin" && "NWCHEM_MODULES" == "tce" ]]; then
    do_largeqas=0
 fi
 case "$ARMCI_NETWORK" in
    MPI-PR)
        nprocs=3
        case "$os" in
            Darwin)
                do_largeqas=0
		case "$MPI_IMPL" in
		    openmpi)
			export MPIRUN_NPOPT="-mca mpi_yield_when_idle 0 --oversubscribe -np "
			;;
		esac
            ;;
        esac
	;;
    MPI-MT)
        do_largeqas=0
        ;;
    MPI-PT)
        do_largeqas=0
        ;;
    MPI3)
        case "$os" in
            Darwin)
                do_largeqas=0
		;;
	esac
	;;
 esac
 if [[ "$NWCHEM_MODULES" == "tce" ]]; then
   cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs tce_n2 tce_ccsd_t_h2o tce_h2o_eomcc
    if  [[ "$do_largeqas" == 1 ]]; then
	cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs tce_ipccsd_f2 tce_eaccsd_ozone
    fi
 else
     cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs $nprocs dft_he2+ prop_mep_gcube
     head -2 $NWCHEM_TOP/QA/testoutputs/dft_he2+.out
     tail -20 $NWCHEM_TOP/QA/testoutputs/dft_he2+.out
     cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs cosmo_h2o_dft  
     cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs pspw 
     if  [[ "$do_largeqas" == 1 ]]; then
       cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs dft_siosi3 h2o_opt
       cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs tddft_h2o h2o2-response
       cd $NWCHEM_TOP/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs pspw_md
     fi
 fi
