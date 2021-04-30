#!/usr/bin/env bash
rm -f dftd3.f nwpwxc_vdw3a.F
if [[ -f "dftd3.tgz" ]]; then
    echo "using existing" dftd3.tgz
else
    echo "downloading"  dftd3.tgz
    wget https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/dftd3.tgz
fi    
tar xzf dftd3.tgz dftd3.f
mv dftd3.f nwpwxc_vdw3a.F
patch -p0 < nwpwxc_vdw3a.patch

