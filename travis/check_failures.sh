#!/bin/bash
tail -20 $TRAVIS_BUILD_DIR/src/tools/build/config.log
tail -10 $TRAVIS_BUILD_DIR/src/tools/build/comex/config.log
grep -A 2 -B 2 -i error $TRAVIS_BUILD_DIR/src/make.log 
head -200 $TRAVIS_BUILD_DIR/src/make.log
tail -200 $TRAVIS_BUILD_DIR/src/make.log
tail -200 $TRAVIS_BUILD_DIR/src/6log
grep -i tce_energy $TRAVIS_BUILD_DIR/src/6log
# tail output files on failures
if [[  "NWCHEM_MODULES" == "tce" ]]; then
    head -2 $TRAVIS_BUILD_DIR/QA/testoutputs/tce_n2.out
    tail -70 $TRAVIS_BUILD_DIR/QA/testoutputs/tce_n2.out
    tail -70 $TRAVIS_BUILD_DIR/QA/testoutputs/tce_h2o_eomcc.out
    tail -70 $TRAVIS_BUILD_DIR/QA/testoutputs/tce_ccsd_t_h2o
else    
    grep d= $TRAVIS_BUILD_DIR/QA/testoutputs/dft_he2+.out
    head -2 $TRAVIS_BUILD_DIR/QA/testoutputs/dft_he2+.out
    tail -40 $TRAVIS_BUILD_DIR/QA/testoutputs/dft_he2+.out
    tail -40 $TRAVIS_BUILD_DIR/QA/testoutputs/prop_mep_gcube.out
    tail -60 $TRAVIS_BUILD_DIR/QA/testoutputs/dft_siosi3.out
    grep @ $TRAVIS_BUILD_DIR/QA/testoutputs/h2o_opt.out
    tail -60 $TRAVIS_BUILD_DIR/QA/testoutputs/h2o_opt.out
    tail -60 $TRAVIS_BUILD_DIR/QA/testoutputs/tddft_h2o.out
    tail -60 $TRAVIS_BUILD_DIR/QA/testoutputs/h2o2-response.out
    tail -490 $TRAVIS_BUILD_DIR/QA/testoutputs/h2o2-response.out
    tail -60 $TRAVIS_BUILD_DIR/QA/testoutputs/pspw.out
    tail -60 $TRAVIS_BUILD_DIR/QA/testoutputs/pspw_md.out
    grep 'Total PSPW energy' $TRAVIS_BUILD_DIR/QA/testoutputs/pspw_md.out
    fi
