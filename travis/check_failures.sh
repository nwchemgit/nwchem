#!/usr/bin/env bash
if [[ -z "$TRAVIS_BUILD_DIR" ]] ; then
    TRAVIS_BUILD_DIR=$(pwd)
fi
source $TRAVIS_BUILD_DIR/travis/nwchem.bashrc
env
head -4000 $TRAVIS_BUILD_DIR/src/make.log
grep -A 2 -B 2 -i error $TRAVIS_BUILD_DIR/src/make.log 
cat $TRAVIS_BUILD_DIR/src/tools/build/config.log
cat $TRAVIS_BUILD_DIR/src/tools/build/comex/config.log
tail -4000 $TRAVIS_BUILD_DIR/src/make.log
ls -lrt $TRAVIS_BUILD_DIR/src/libext/lib/
echo '###### OpenBLAS make.log ####'
grep -i  gemm_incopy $TRAVIS_BUILD_DIR/src/make.log
grep  NO_LAPACKE $TRAVIS_BUILD_DIR/src/make.log 
echo '###### end of make.log ####'
if [[ "$USE_64TO32" == "y" ]]; then
    tail -200 $TRAVIS_BUILD_DIR/src/6log
    grep -i tce_energy $TRAVIS_BUILD_DIR/src/6log
fi
# tail output files on failures
check_file () {
    file=$TRAVIS_BUILD_DIR/QA/testoutputs/$1.out
    if [ -f $file ] ; then
	echo ============================================================
	echo $file
	head -2 $file
	echo ============================================================
	tail -260 $file
	echo ============================================================
	grep -s 'l DFT energy' $file |tail
	echo ============================================================
	grep -s 'Total PSPW energy' $file |tail
    fi
}
    
check_file tce_n2
check_file tce_h2o_eomcc
check_file tce_ccsd_t_h2o
check_file dft_he2+
check_file prop_mem_gcube
check_file dft_siosi3
check_file h2o_opt
check_file tddft_h2o
check_file h2o2-response
check_file pspw
check_file pspw_md
check_file aump2
# stuff for OpenBLAS issues on macos
echo 'looking for cgemm_incopy.o'
find $TRAVIS_BUILD_DIR/src/libext -name cgemm_incopy.o
echo ' looking for _cgemm_incopy symbol in cgemm_incopy.o'
nm `find $TRAVIS_BUILD_DIR/src/libext -name cgemm_incopy.o` 
echo ' done OpenBLAS debugging '
