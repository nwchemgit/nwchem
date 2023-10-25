#!/usr/bin/env bash
sudo apt-get install -y mpich libmpich-dev
export NWCHEM_TOP=`pwd`
export USE_MPI=1
cd src
make nwchem_config NWCHEM_MODULES="all python"
make  64_to_32 CONVERT_ALL=y
make  32_to_64
rm -f diff.out
git diff -U0 . >& diff.out
if [ $(wc -l diff.out | cut -d " " -f 1) != 0 ]; then
    cat diff.out
    echo "********** check_64_to_32 *********"
    echo "********** found missing files ****"
    echo "********** from USES_BLAS *********"
    grep 'diff --git' diff.out  | cut -d ' ' -f 4 | sed -e "s/b\/src/src/" 
    echo "***********************************"
    exit 1
else
    echo "********** check_64_to_32 *********"
    echo "********** found no missing files ****"
    echo "********** from USES_BLAS *********"
fi
