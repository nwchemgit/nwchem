#!/bin/bash
#
# This script works out the revision number of the GA source
# code. It writes the resulting data in suboutine that can be used
# to query this information. As NWChem can build with different 
# versions of GA the script takes one command line argument to
# figure out of which GA directory it has to establish the revision
# number.
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
my_svnversion=`which svn`
ga_dir="$1"
cd "$path"
if [ -f "${my_svnversion}" ] ; then
  # svnversion exists, but does .svn?
  if [ -d "../tools/${ga_dir}/.svn" ] ; then
    # .svn exists too
    revision=`${my_svnversion} info "../tools/${ga_dir}" | grep Revision:`
    revision=`echo ${revision} | sed 's/Revision: //'`
    echo "      subroutine util_ga_version(garev)" > util_ga_version.F
    echo "      implicit none" >> util_ga_version.F
    echo "      character*(*) garev" >> util_ga_version.F
    echo "      garev=\"${revision}\"" >> util_ga_version.F
    echo "      end" >> util_ga_version.F
  else
    if [ ! -f util_ga_version.F ] ; then
      echo "      subroutine util_ga_version(garev)" > util_ga_version.F
      echo "      implicit none" >> util_ga_version.F
      echo "      character*(*) garev" >> util_ga_version.F
      echo "      garev=\"N/A\"" >> util_ga_version.F
      echo "      end" >> util_ga_version.F
    fi
  fi
else
  if [ ! -f util_ga_version.F ] ; then
    echo "      subroutine util_ga_version(garev)" > util_ga_version.F
    echo "      implicit none" >> util_ga_version.F
    echo "      character*(*) garev" >> util_ga_version.F
    echo "      garev=\"N/A\"" >> util_ga_version.F
    echo "      end" >> util_ga_version.F
  fi
fi
