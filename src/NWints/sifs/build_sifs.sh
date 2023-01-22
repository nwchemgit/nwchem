#!/usr/bin/env bash
TGZ=columbus-master.tar.gz
[ -f "$TGZ" ] && gunzip -t "$TGZ" > /dev/null && USE_TGZ=1
if [[ "$USE_TGZ" == 1 ]]; then
    echo "using existing"  "$TGZ"
else
rm -rf columbu* *.F sifs.patched
tries=0 ; until [ "$tries" -ge 10 ] ; do \
curl > "$TGZ" https://gitlab.com/api/v4/projects/36816383/repository/archive?path=Columbus/source/gcfci/colib/sifs \
&& curl > bummer.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/humanio/bummer.F?inline=false\
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
