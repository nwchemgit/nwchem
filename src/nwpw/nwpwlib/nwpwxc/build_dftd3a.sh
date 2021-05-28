#!/usr/bin/env bash
rm -f dftd3.f nwpwxc_vdw3a.F
URL1="https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/"
URL2="https://web.archive.org/web/20210527062154if_/https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/"
declare -a urls=("$URL1" "$URL1" "$URL1" "$URL2" "$URL2")
TGZ=dftd3.tgz
if [[ -f "$TGZ" ]]; then
    echo "using existing" "$TGZ"
else
    echo "downloading"  "$TGZ"
    CURL_YES=`curl  -O 2>&1 | head -1  | awk ' /URL/ {print "Y";exit};{print "N"}'`
    tries=1 ; until [ "$tries" -ge 5 ] ; do
    if [ $CURL_YES = "Y1" ]; then
	curl -f -L ${urls[$tries]}/"$TGZ" -o "$TGZ" && break
    else
	wget ${urls[$tries]}/"$TGZ" && break
    fi
    tries=$((tries+1)) ; echo attempt no.  $tries    ; sleep 5 ;  done
fi
if [[ ! -f "$TGZ" ]]; then
    echo "download failed"
    exit 1
fi
tar xzf dftd3.tgz dftd3.f
mv dftd3.f nwpwxc_vdw3a.F
patch -p0 < nwpwxc_vdw3a.patch

