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
cd ../..
rm -rf cpuinfo
