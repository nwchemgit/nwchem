#!/usr/bin/env bash
TGZ=columbus-master.tar.gz
[ -f "$TGZ" ] && gunzip -t "$TGZ" > /dev/null && USE_TGZ=1
USE_TGZ=0
if [[ "$USE_TGZ" == 1 ]]; then
    echo "using existing"  "$TGZ"
else
  rm -rf columbu* $(ls *.F |sed 's/sifs_stubs.F//')  sifs.patched
tries=0 ; until [ "$tries" -ge 10 ] ; do \
curl > "$TGZ" https://gitlab.com/api/v4/projects/36816383/repository/archive?path=Columbus/source/gcfci/colib/sifs \
&& curl > bummer.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/humanio/bummer.F?inline=false\
&& curl > bummer.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/humanio/bummer.F?inline=false\
&& curl > aiopen.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/basicio/aiopen.F?inline=false\
&& curl > airead.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/basicio/airead.F?inline=false\
&& curl > aiwait.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/basicio/aiwait.F?inline=false\
&& curl > aiwrit.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/basicio/aiwrit.F?inline=false\
&& curl > seqwbf.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/basicio/seqwbf.F?inline=false\
&& curl > seqrbf.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/basicio/seqrbf.F?inline=false\
&& curl > plab1.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/pack/plab1.F?inline=false\
&& curl > plab8.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/pack/plab8.F?inline=false\
&& curl > plab16.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/pack/plab16.F?inline=false\
&& curl > plab32.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/pack/plab32.F?inline=false\
&& curl > ulab1.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/pack/ulab1.F?inline=false\
&& curl > ulab8.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/pack/ulab8.F?inline=false\
&& curl > ulab16.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/pack/ulab16.F?inline=false\
&& curl > ulab32.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/pack/ulab32.F?inline=false\
&& curl > izero.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/linear_algebra/izero.F?inline=false\
	    && break ;\
					 tries=$((tries+1)) ; echo attempt no.  $tries    ; sleep 2 ;  done
fi
[ -f "$TGZ" ] && gunzip -t "$TGZ" > /dev/null && USE_TGZ=1
if [[ "$USE_TGZ" != 1 ]]; then
    echo "download failed"
    echo 
    echo "if internet connectivity is missing"
    echo "set NO_SIFS=1"
    echo 
    exit 1
fi
tar xzf "$TGZ"
mv $(find columbus-master*sifs -name "*F") .
rm -rf columbus-master*sifs
if [[ -f "sifs.patched" ]]; then
    echo 'source already patched'
else
rm -f ibummr.patch
cat > ibummr.patch <<EOF
--- bummer.F.org	2022-12-01 11:34:21.510551238 -0800
+++ bummer.F	2022-12-01 11:32:28.997830864 -0800
@@ -201,7 +201,7 @@
 c
 c     # initialization...
 c
-      entry ibummr( iunit )
+c      entry ibummr( iunit )
 c
 c     # save the listing unit for use later.
 c
EOF
patch -p0 -s -N < ibummr.patch
rm -f ibummr.patch
rm -f siftdy.patch
cat > siftdy.patch <<EOF
--- siftdy.F.org	2022-12-02 14:50:03.807903254 -0800
+++ siftdy.F	2022-12-02 14:54:13.969471516 -0800
@@ -282,7 +282,7 @@
 c
 c     # default code: just return a dummy string.
 c
-      chrtdy = site // 'Machine=?  ??:??:?? ??-???-??'
+      chrtdy = site // 'Machine=unknown'
 #endif 
 c
       return
EOF
patch -p0 -s -N < siftdy.patch
rm -f siftdy.patch
echo yes > sifs.patched
fi

# dcopy is cursed, won't work with COLUMBUS
sed -i -e 's/dcopy/cdcopy/g' ./sif[d,e][1,2].F

# SO effective density matrix has btype 35.  Needed until this change is made in COLUMBUS
sed -i -e 's/ibtypmx=34/ibtypmx=40/g' ./sifr1n.F
sed -i -e 's/btypmx=20/btypmx=40/g' ./sifr1x.F
