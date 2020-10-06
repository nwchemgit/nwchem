#!/bin/bash 
# Exit on error
set -ev
# source env. variables
source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
# check if nwchem binary has been cached
if [[ -f "$NWCHEM_EXECUTABLE" ]] ; then
    EXTRA_BUILD=0
else
    echo 'Cached NWChem binary not found, recompiling'
    $TRAVIS_BUILD_DIR/travis/config_nwchem.sh 
    $TRAVIS_BUILD_DIR/travis/compile_nwchem.sh
    EXTRA_BUILD=1
fi
os=`uname`
arch=`uname -m`
export NWCHEM_BASIS_LIBRARY=$TRAVIS_BUILD_DIR/.cachedir/files/libraries/
export NWCHEM_NWPW_LIBRARY=$TRAVIS_BUILD_DIR/.cachedir/files/libraryps/
 nprocs=2
 do_largeqas=1

 if [[ "$EXTRA_BUILD" == "1" ]]; then
     do_largeqas=0
 fi

 if [[ "$os" == "Linux" && "$MPI_IMPL" == "mpich" ]]; then
    export MPIRUN_PATH=/usr/bin/mpirun.mpich
 fi
 if [[ "$os" == "Darwin" && "NWCHEM_MODULES" == "tce" ]]; then
    do_largeqas=0
 fi
 case "$ARMCI_NETWORK" in
    MPI-PR)
	nprocs=$(( nprocs + 1 ))
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
    SOCKETS)
        do_largeqas=0
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
 echo === ls binaries cache ===
 ls -lrt $TRAVIS_BUILD_DIR/.cachedir/binaries/$NWCHEM_TARGET/ || true
 echo =========================
 if [[ "$NWCHEM_MODULES" == "tce" ]]; then
   cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs tce_n2 tce_ccsd_t_h2o tce_h2o_eomcc
   cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs ducc_be
 if  [[ "$do_largeqas" == 1 ]]; then
	cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs tce_ipccsd_f2 tce_eaccsd_ozone
    fi
 else
     cd $TRAVIS_BUILD_DIR/QA && ./runtests.mpi.unix procs $nprocs dft_he2+ prop_mep_gcube
     cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs cosmo_h2o_dft  
     if [[ "$USE_SIMINT" != "1" ]] ; then
	cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs pspw
     fi
     if [[ "$NWCHEM_MODULES" == "tinyqmpw python" ]]; then
       cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs pyqa3
     fi
     if  [[ "$do_largeqas" == 1 ]]; then
       cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs dft_siosi3 h2o_opt
       cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs tddft_h2o h2o2-response
       cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs dft_scan
       cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs dft_ncap
       cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs dft_ch3_h2o_revm06
       cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs dft_smear
       cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs n2_ccsd h2mp2 auh2o aump2
       cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs pspw_md
       cd $TRAVIS_BUILD_DIR/QA && USE_SLEEPLOOP=1 ./runtests.mpi.unix procs $nprocs ch3radical_unrot
     fi
 fi
