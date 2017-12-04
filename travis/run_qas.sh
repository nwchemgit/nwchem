#!/bin/bash -f
# source env. variables
 source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
 if [[ "$NWCHEM_MODULES" == "tce" ]]; then
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs 2 tce_n2 tce_ccsd_t_h2o tce_h2o_eomcc
   echo "diff output for tce_h2o_eomcc"
   diff $NWCHEM_TOP/QA/testoutputs/tce_h2o_eomcc.ok.out.nwparse $NWCHEM_TOP/QA/testoutputs/tce_h2o_eomcc.out.nwparse
   ls -l $NWCHEM_TOP/QA/testoutputs/tce_h2o_eomcc.out.nwparse $NWCHEM_TOP/QA/testoutputs/tce_h2o_eomcc.out.nwparse \
   $NWCHEM_TOP/QA/testoutputs/tce_h2o_eomcc.oute $NWCHEM_TOP/QA/testoutputs/tce_h2o_eomcc.out
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs 2 tce_ipccsd_f2
   echo "diff output for tce_ipccsd_f2"
   diff $NWCHEM_TOP/QA/testoutputs/tce_ipccsd_f2.ok.out.nwparse $NWCHEM_TOP/QA/testoutputs/tce_ipccsd_f2.out.nwparse
 else
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs 2 dft_siosi3
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs  2 h2o_opt dft_he2+ cosmo_h2o_dft tddft_h2o h2o2-response
   cd $NWCHEM_TOP/QA && ./runtests.mpi.unix procs  2 pspw 
 fi
