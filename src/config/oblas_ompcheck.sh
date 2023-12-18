#!/usr/bin/env bash
# oblas_dir directory where libopenblas is located
# oblas_name "openblas" suffix in lib"openblas".
# e.g. for debian pkg libopenblas64-pthread-dev
# oblas_dir=/usr/lib/x86_64-linux-gnu/openblas64-pthread/
# oblas_name=openblas64
#echo "USE_OPENMP is equal to " $USE_OPENMP
mylog=/tmp/mylog.txt
NWCHEM_TOP=$1
file_check=$NWCHEM_TOP/src/oblas_ompcheck_done.txt
if [ -f $file_check ]; then
#    echo $file_check present >> $mylog
    cat $file_check
    exit 0
fi
echo BLASOPT is $BLASOPT >> $mylog
if [[  -z "${BLASOPT}" ]]; then
    echo BLASOPT is empty >> $mylog
    echo 0 >  $NWCHEM_TOP/src/oblas_ompcheck_done.txt
    exit 0
else
    echo BLASOPT is not empty @$BLASOPT@ >> $mylog
fi
#extract oblas_dir & oblas_name from BLASOPT
# 
oblas_dir=$(echo $BLASOPT | awk 'sub(/.*-L */,""){f=1} f{if ( sub(/ * .*/,"") ) f=0; print}')
oblas_name=$(echo $BLASOPT | awk 'sub(/.* -l */,""){f=1} f{if ( sub(/ * .*/,"") ) f=0; print}')
echo oblas_dir is $oblas_dir  oblas_name is $oblas_name >> $mylog
find $oblas_dir -name "lib*openblas*.*"|| true
if [ $(uname -s) == 'Darwin' ]; then
    MYLDD='otool -L'
    SOSUFFIX=dylib
else
    MYLDD=ldd
    SOSUFFIX=so
fi
# check first against clang libomp
if [ -f "$oblas_dir/lib$oblas_name.$SOSUFFIX" ]; 
then
    gotomp=$($MYLDD $oblas_dir/lib$oblas_name.$SOSUFFIX | grep libomp | wc -l )
# next check against gcc libgomp
    if [ $gotomp -eq 0 ]
    then
	gotomp=$($MYLDD $oblas_dir/lib$oblas_name.$SOSUFFIX | grep libgomp | wc -l )
    fi
else
    gotomp=0
fi
echo gotomp $gotomp >> $mylog
#conda packages might use OpenMP to thread OpenBLAS
if [ $gotomp -ne 0 ]
then
    echo openblas built with OpenMP >> $mylog
    export OPENBLAS_USES_OPENMP=1
    unset USE_OPENMP
else
    echo openblas built without OpenMP >> $mylog
fi
echo $gotomp >  $NWCHEM_TOP/src/oblas_ompcheck_done.txt
echo $gotomp 
exit 0
