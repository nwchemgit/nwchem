#!/usr/bin/env bash
if [ ! -f ${NWCHEM_TOP}/src/mpi_include.txt ]; then
    ${NWCHEM_TOP}/src/tools/guess-mpidefs --mpi_include > ${NWCHEM_TOP}/src/mpi_include.txt
fi
cat  ${NWCHEM_TOP}/src/mpi_include.txt 
