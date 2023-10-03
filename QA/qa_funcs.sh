#!/usr/bin/env bash
get_nwchem_top(){
    if [[ -z "${NWCHEM_TOP}" ]]; then
	DIRQA=`dirname "$0"`
	MYPWD=`pwd`
	NWCHEM_TOP=`echo ${MYPWD}/${DIRQA} | sed -e 's/\/QA.*//' `
    fi
    echo $NWCHEM_TOP
}

get_nwchem_target(){
    if [[ -z "${NWCHEM_TARGET}" ]]; then
	UNAME_S=$(uname -s)
	if [[ ${UNAME_S} == Linux ]]; then
	    NWCHEM_TARGET=LINUX64
	elif [[ ${UNAME_S} == Darwin ]]; then
	    NWCHEM_TARGET=MACX64
	else
	    echo
	    echo You must define NWCHEM_TARGET in your environment to be the name
	    echo of the machine you wish to build for ... for example
	    echo     export NWCHEM_TARGET=SOLARIS
	    echo Known targets are SOLARIS, ...
	    echo See the INSTALL instructions for a complete list
	    echo ${UNAME_S}
	    exit 1
	fi
    fi
    echo $NWCHEM_TARGET
}    
get_nwchem_executable(){
    NWCHEM_TOP=$(get_nwchem_top)
    NWCHEM_TARGET=$(get_nwchem_target)
    if [[ $(echo $NWCHEM_EXECUTABLE|cut -c 1-11) == 'singularity' ]] || [[ $(echo $NWCHEM_EXECUTABLE|cut -c 1-9) == 'apptainer' ]] || [[ -f $NWCHEM_EXECUTABLE ]]; then
	NWCHEM=$NWCHEM_EXECUTABLE
    else
	NWCHEM=${NWCHEM_TOP}/bin/${NWCHEM_TARGET}/nwchem
    fi
    echo $NWCHEM
}

