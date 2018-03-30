#!/usr/bin/env bash
if [ -f "$0" ] ; then
   # The first item on the command line is an actual file so the 
   # script must have been specified including the path.
   path="`dirname \"$0\"`"
else
   # The first item on the command line is not a file so script
   # it must have been found in PATH.
   path="`which \"$0\"`"
   path="`dirname \"$path\"`"
fi
if [[ -z "${EXTERNAL_GA_PATH}" ]]; then
    GA_PATH=${NWCHEM_TOP}/src/tools/install
else
    GA_PATH=${EXTERNAL_GA_PATH}
fi
gotscalapack=`${GA_PATH}/bin/ga-config  --use_scalapack |  awk ' /1/  {print ".true.";exit};{print ".false."}'`

cd "$path"
rm -f util_scalapack_info.F
echo ggog $gotscalapack
    echo "      logical function util_scalapack_info()" > util_scalapack_info.F
    echo "      implicit none" >> util_scalapack_info.F
    echo "      util_scalapack_info=${gotscalapack}" >> util_scalapack_info.F
    echo "      end" >> util_scalapack_info.F
