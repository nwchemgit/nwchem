#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs 2 dft_siosi3
 cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs  2 h2o_opt dft_he2+ cosmo_h2o_dft tddft_h2o h2o2-response
