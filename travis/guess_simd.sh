#!/usr/bin/env bash
UNAME_S=$(uname -s)
arch=`uname -m`
if [[ ${UNAME_S} == Linux ]]; then
    CPU_FLAGS=$(lscpu | grep -i flags | sort | uniq)
elif [[ ${UNAME_S} == Darwin ]]; then
    CPU_FLAGS=$(/usr/sbin/sysctl -n machdep.cpu.features)
    if [[ "$arch" == "x86_64" ]]; then
	CPU_FLAGS_2=$(/usr/sbin/sysctl -n machdep.cpu.leaf7_features)
    fi
fi
if [[ $(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /clzero/{print "Y"}') == "Y" ]]; then
    echo "zen"
    exit 0
elif [[ $(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /avx512f/{print "Y"}') == "Y" ]]; then
    echo "avx512"
    exit 0
elif [[ $(echo ${CPU_FLAGS} | tr  'A-Z' 'a-z'| awk ' /avx2/   {print "Y"}') == "Y" ]]; then
    echo "avx2"
    exit 0
elif [[ $(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /avx/    {print "Y"}') == "Y" ]]; then
    echo "avx"
    exit 0
elif [[ $(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /sse2/   {print "Y"}') == "Y" ]]; then
    echo "sse2"
    exit 0
elif [[ $arch == "arm64" || $arch == "aarch64" ]]; then
    if [[ $(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /sve2/   {print "Y"}') == "Y" ]]; then
	echo "sve2"
	exit 0
    elif [[ $(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /sve/   {print "Y"}') == "Y" ]]; then
	echo "sve"
	exit 0
    elif [[ $(echo ${CPU_FLAGS}   | tr  'A-Z' 'a-z'| awk ' /asimddp/   {print "Y"}') == "Y" ]]; then
	echo "asimddp"
	exit 0
    else
	echo "arm64"
	exit 0
    fi
else
    echo "unknown"
fi
