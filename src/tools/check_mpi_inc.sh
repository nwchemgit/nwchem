#!/usr/bin/env bash
echo `date` >> /tmp/edo/log
if [ ! -f ${NWCHEM_TOP}/src/mpi_include.txt ]; then
    ${NWCHEM_TOP}/src/tools/guess-mpidefs --mpi_include > ${NWCHEM_TOP}/src/mpi_include.txt
    echo created file  >> /tmp/edo/log
fi
/usr/bin/cat  ${NWCHEM_TOP}/src/mpi_include.txt >> /tmp/edo/log
/usr/bin/cat  ${NWCHEM_TOP}/src/mpi_include.txt 
