#!/usr/bin/env bash
# convert all source files
find ${NWCHEM_TOP}/src \
     -type d \( -path tools -o -path config  \
  -o -path lapack  -o -path blas \
     -o -path libext -o -name libsimint_source \) -prune  \
-o -type f \( -iname "*.f"  -o  -iname "*.c"  \) \
  | xargs -n 50 ${NWCHEM_TOP}/src/config/64_to_32 
