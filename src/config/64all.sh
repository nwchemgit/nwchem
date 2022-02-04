#!/usr/bin/env bash
# convert all source files
find ${NWCHEM_TOP}/src/NWints \
	 -type d \( -name "blas*" \
	 -o -name "lapack*" \
	 -o -name "config*" \
	 -o -name "libext*" \
	 -o -name "tools*" \
	 -o -name "lapack*" \
	 -o -name "libsimint_source*" \
	 -o -name "64to32blas*" \) -prune \
-o -type f \( -iname "*.f"  -o  -iname "*.c"  \) \
  | xargs -n 50 ${NWCHEM_TOP}/src/config/64_to_32
