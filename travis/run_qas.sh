#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 if [[ "$NWCHEM_MODULES" == "tce" ]]; then
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs 2 tce_n2
 else
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs 2 dft_siosi3
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs  2 h2o_opt dft_he2+ cosmo_h2o_dft tddft_h2o h2o2-response
 fi
