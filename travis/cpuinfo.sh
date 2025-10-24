#!/usr/bin/env bash
git clone -q https://github.com/pytorch/cpuinfo
cd cpuinfo
mkdir build && cd build
CMAKE_POLICY_VERSION_MINIMUM=3.5 cmake  .. >& cmake.log
tail -20 cmake.log
make -j3 >& make.log
tail -3 make.log
./isa-info
./cpu-info
./cpu-info |grep -A 1 Microa|tail -1 |cut -d ' ' -f 2-3 |  sed -e "s/ //"
rm -f /tmp/microarch_$(id -u).txt
./cpu-info |grep -A 1 Microa|tail -1 |cut -d ' ' -f 2-3 |  sed -e "s/ //"> /tmp/microarch_$(id -u).txt
cd ../..
if [[ $(uname -s) == Linux ]]; then /usr/bin/lscpu ; fi
