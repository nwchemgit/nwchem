#!/usr/bin/env bash
git clone https://github.com/pytorch/cpuinfo
cd cpuinfo
mkdir build && cd build
cmake ..
make -j3
./isa-info
./cpu-info
cd ../..
rm -rf cpuinfo
