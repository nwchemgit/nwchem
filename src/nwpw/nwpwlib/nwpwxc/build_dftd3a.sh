#!/usr/bin/env bash
rm -f dftd3.f nwpwxc_vdw3a.F
URL="https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/"
TGZ=dftd3.tgz
if [[ -f "$TGZ" ]]; then
    echo "using existing" "$TGZ"
else
    echo "downloading"  "$TGZ"
    CURL_YES=`curl  -O 2>&1 | head -1  | awk ' /URL/ {print "Y";exit};{print "N"}'`
    if [ $CURL_YES = "Y" ];	then
	curl -L "$URL"/"$TGZ" -o "$TGZ"
    else
	wget "$URL"/"$TGZ"
    fi
fi    
tar xzf dftd3.tgz dftd3.f
mv dftd3.f nwpwxc_vdw3a.F
patch -p0 < nwpwxc_vdw3a.patch

