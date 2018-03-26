#!/usr/bin/env bash
#
# This script works out the revision number of the NWChem source
# code. It writes the resulting data in suboutine that can be used
# to query this information.
#
# We need 2 things for this operation to complete successfully:
# 1. svn needs to be available on this machine
# 2. the .svn directories need to be present in this source code
# If both these requirements are satisfied we will always overwrite
# the revision information with a current version.
# If either of these requirements is not satisfied we still need to
# make sure there is a valid version subroutine to ensure the code
# will compile. If such a routine already exists we do nothing as 
# we lack the tools to do better than whatever is in that routine
# already. If such a routine does not exist we create one but
# setting the version number to blank (this is the best we can do).
#
# First find out where this script actually lives so we can create
# the appropriate Fortran file in the right location.
#
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
my_gitversion=`which git`
cd "$path"
if [ -f "${my_gitversion}" ] ; then
  # gitversion exists, but is the code under git?
WCBRANCH=`${my_gitversion}  describe  --always   | wc -l`
if [ ${WCBRANCH} -ne 0 ]; then
    revision=`${my_gitversion}  describe  --always`
    echo "      subroutine util_nwchem_version(nwrev)" > util_nwchem_version.F
    echo "      implicit none" >> util_nwchem_version.F
    echo "      character*(*) nwrev" >> util_nwchem_version.F
    echo "      nwrev=\"${revision}\"" >> util_nwchem_version.F
    echo "      end" >> util_nwchem_version.F
  else
    if [ ! -f util_nwchem_version.F ] ; then
      echo "      subroutine util_nwchem_version(nwrev)" > util_nwchem_version.F
      echo "      implicit none" >> util_nwchem_version.F
      echo "      character*(*) nwrev" >> util_nwchem_version.F
      echo "      nwrev=\"N/A\"" >> util_nwchem_version.F
      echo "      end" >> util_nwchem_version.F
    fi
  fi
else
  if [ ! -f util_nwchem_version.F ] ; then
    echo "      subroutine util_nwchem_version(nwrev)" > util_nwchem_version.F
    echo "      implicit none" >> util_nwchem_version.F
    echo "      character*(*) nwrev" >> util_nwchem_version.F
    echo "      nwrev=\"N/A\"" >> util_nwchem_version.F
    echo "      end" >> util_nwchem_version.F
  fi
fi
