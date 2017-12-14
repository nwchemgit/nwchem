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
my_gitversion=`which git`
if [ $# -eq 0 ]
  then
# if no arg supplied,
# try to guess from tools/build
    if [  -f $NWCHEM_TOP/src/tools/build/config.log ] ; then

   ga_dir=`head  -7 $NWCHEM_TOP/src/tools/build/co*log|tail -1 |cut -d '/' -f2`
   else
    echo "ga_dir argument not supplied, will write N/A revision"
   fi
    echo $ga_dir
else
    ga_dir="$1"
fi

if [ -d "$NWCHEM_TOP/src/tools/${ga_dir}/.git" ] ; then
    cd $NWCHEM_TOP/src/tools/${ga_dir}
    if [ -f "${my_gitversion}" ] ; then
	# gitversion exists, but is the ga_dir under git?
	revision="N/A"
	GITBRANCH=`${my_gitversion} describe --tags 2> /dev/null| wc -l`
	if [ ${GITBRANCH} -ne 0 ]; then
	    # 
	    revision=`${my_gitversion} describe --tags`
	fi
    fi
else
    revision=${ga_dir}
fi
cd $NWCHEM_TOP/src/util
    echo "      subroutine util_ga_version(garev)" > util_ga_version.F
    echo "      implicit none" >> util_ga_version.F
    echo "      character*(*) garev" >> util_ga_version.F
    echo "      garev=\"${revision}\"" >> util_ga_version.F
    echo "      end" >> util_ga_version.F
