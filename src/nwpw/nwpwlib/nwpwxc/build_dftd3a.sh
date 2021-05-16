#!/usr/bin/env bash
rm -f dftd3.f nwpwxc_vdw3a.F
URL="https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/"
TGZ=dftd3.tgz
if [[ -f "$TGZ" ]]; then
    echo "using existing" "$TGZ"
else
    echo "downloading"  "$TGZ"
    CURL_YES=`curl  -O 2>&1 | head -1  | awk ' /URL/ {print "Y";exit};{print "N"}'`
    tries=0 ; until [ "$tries" -ge 5 ] ; do
    if [ $CURL_YES = "Y" ]; then
	curl -L "$URL"/"$TGZ" -o "$TGZ" && break
    else
	wget "$URL"/"$TGZ" && break
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

