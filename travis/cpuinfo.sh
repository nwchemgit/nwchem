#!/usr/bin/env bash
git clone -q https://github.com/pytorch/cpuinfo
cd cpuinfo
mkdir build && cd build
cmake .. >& cmake.log
tail -3 cmake.log
make -j3 >& make.log
tail -3 make.log
./isa-info
./cpu-info
cd ../..
rm -rf cpuinfo
