#!/usr/bin/env bash
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

if [[ -z "${NWCHEM_TOP}" ]]; then
    DIRUTIL=`dirname "$0"`
    MYPWD=`pwd`
    NWCHEM_TOP=`echo ${MYPWD}/${DIRUTIL} | sed -e 's/\/src.*//' `
fi
# check if EXTERNAL_GA_PATH is set

    revision="N/A"
if [[ -z "${EXTERNAL_GA_PATH}" ]]; then
    # try to guess from tools/build
    if [  -f $NWCHEM_TOP/src/tools/build/config.log ] ; then

	ga_dir=`head  -7 $NWCHEM_TOP/src/tools/build/co*log|tail -1 |cut -d '/' -f2`
	if [ -d "$NWCHEM_TOP/src/tools/${ga_dir}/.git" ] ; then
            my_gitversion=`which git`
	    cd $NWCHEM_TOP/src/tools/${ga_dir}
	    if [ -f "${my_gitversion}" ] ; then
		# gitversion exists, but is the ga_dir under git?
		GITBRANCH=`${my_gitversion} describe --tags 2> /dev/null| wc -l`
		if [ ${GITBRANCH} -ne 0 ]; then
		    # 
		    revision=`${my_gitversion} describe --tags`
		fi
	    fi
	fi
    fi
else
    revision=`${EXTERNAL_GA_PATH}/bin/ga-config --version`
fi

if [ "$revision" == "N/A" ] ; then
#no .git information
    if [ -f ${NWCHEM_TOP}/src/tools/install/bin/ga-config ] ; then
	revision=`${NWCHEM_TOP}/src/tools/install/bin/ga-config --version`
    fi
fi

cd $NWCHEM_TOP/src/util
    echo "      subroutine util_ga_version(garev)" > util_ga_version.F
    echo "      implicit none" >> util_ga_version.F
    echo "      character*(*) garev" >> util_ga_version.F
    echo "      garev=\"${revision}\"" >> util_ga_version.F
    echo "      end" >> util_ga_version.F
