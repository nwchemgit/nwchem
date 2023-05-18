#!/usr/bin/env bash
TGZ=columbus-master.tar.gz
[ -f "$TGZ" ] && gunzip -t "$TGZ" > /dev/null && USE_TGZ=1
USE_TGZ=0
if [[ "$USE_TGZ" == 1 ]]; then
    echo "using existing"  "$TGZ"
else
rm -rf columbu* *.F sifs.patched
tries=0 ; until [ "$tries" -ge 10 ] ; do \
curl > "$TGZ" https://gitlab.com/api/v4/projects/36816383/repository/archive?path=Columbus/source/colib \
&& curl > izero.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/linear_algebra/izero.F?inline=false\
&& curl > blaswrapper.f https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/blaswrapper/blaswrapper.f?inline=false\
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
#mv $(find columbus-master*sifs -name "*F") .
#rm -rf columbus-master*sifs
mv $(find columbus-master*colib -name "colib*.f") .
mv $(find columbus-master*colib -wholename "*/colib/sif_data.h") .
#mv columbus-master*colib/sif_data.h .
rm -rf columbus-master*colib
rm -rf newsif
if [[ -f "sifs.patched" ]]; then
    echo 'source already patched'
else
rm -f ibummr.patch
cat > ibummr.patch <<EOF
--- colib7.f.org	2022-12-01 11:34:21.510551238 -0800
+++ colib7.f	2022-12-01 11:32:28.997830864 -0800
@@ -205,7 +205,7 @@
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
cat > siftdy.patch << EOF
--- colib8.f.org
+++ colib8.f
@@ -3160,8 +3160,8 @@ c
 c
 *@ifdef f90 .and. unix
 c     # f90 code with hostnm()
-      integer  hostnm
-      external hostnm
+!      integer  hostnm
+!      external hostnm
       integer :: time(8), ierr
       intrinsic date_and_time, len
       character(len=3), parameter :: month(12) = (/
@@ -3172,11 +3172,11 @@ c     # use the 2-digit year to make more space for the hostname.
       character(len=*), parameter :: cfmt=
      & "(2(i2.2,':'),i2.2,'.',i3.3,i3.2,'-',a3,'-',i2.2)"
 c
-      ierr = hostnm( chrtdy(1:17) )
-      if ( ierr .ne. 0 ) then
+!      ierr = hostnm( chrtdy(1:17) )
+!      if ( ierr .ne. 0 ) then
 c        # pack the error code into chrtdy(:)
-         write( chrtdy(1:17),'(a9,i8.8)' ) 'hostnm()=', ierr
-      endif
+!         write( chrtdy(1:17),'(a9,i8.8)' ) 'hostnm()=', ierr
+!      endif
       chrtdy(18:18) = ' '
       call date_and_time( values=time )
       write( chrtdy(19:40), fmt=cfmt)
EOF
patch -p0 -s -N < siftdy.patch
rm -f siftdy.patch
echo yes > sifs.patched
fi

# dcopy is cursed, won't work with COLUMBUS
sed -i -e 's/dcopy_wr/cdcopy/g' ./colib8.f

# SO effective density matrix has btype 35.  Needed until this change is made in COLUMBUS
sed -i -e 's/btypmx=34/btypmx=40/g' ./colib8.f
sed -i -e 's/btypmx=35/btypmx=40/g' ./colib8.f

#These already exist in NWchem; need to rename to avoid double definition
sed -i -e 's/subroutine\ mxma/subroutine\ mxma_x/g' ./colib1.f
sed -i -e 's/function\ fstrlen/function\ fstrlen_x/g' ./colib7.f
sed -i -e 's/subroutine\ timer/subroutine\ timer_x/g' ./colib3.f
sed -i -e 's/subroutine\ \ icopy(/function\ icopy_x(/g' ./blaswrapper.f
