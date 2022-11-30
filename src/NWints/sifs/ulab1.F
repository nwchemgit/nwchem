      subroutine ulab1( p, u, nuw )
c
c  unpack a bit vector.
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+63)/64)*64 ) are referenced.
c  nuw  = number of unpacked array elements.
c
c  written by ron shepard.
c  version: 5-jul-89
c
      implicit none
      integer nuw
#if defined(CRAY)
      real*8 p(*)
      integer u(nuw)
c
      integer nuw64
c
      nuw64=((nuw+63)/64)*64
      if ( nuw64 .ne. 0 ) call unpack(p,1,u,nuw64)

cevb Need to add definitions for mpack 
cevb for other platforms.  The bitwise functions are
cevb F90/F95 instrinsics

#elif defined(INT64)
        integer p(*), mpack,nuw64,u(64,*),j,l,i,k,plocal ,m2
        parameter (m2=2)
#  if defined(T3E64) || defined(T3D)
         mpack(l,k)= and(shiftr(l,k),1)
#  elif defined(SGIPOWER) || defined(DECALPHA)
         mpack(l,k)= iand(ishft(l,-k),1)
#  elif defined(F90) || defined(F95)
         mpack(l,k)= iand(ishft(l,-k),1)
#  endif 
        nuw64= ((nuw+63)/64)
        call izero(nuw,u,1)
        do 10 i=1, nuw64
         do 11 j=0,63
          u(j+1,i) =mpack(p(i),63-j)
   11    continue
   10    continue
#elif defined(VAX)
c
c  the following code may be used for 32-bit integer, little-endian
c  machines that use vax-type bit operators.
c  30-jul-90  written by ron shepard.
c
      integer p(2,*)
      integer u(*)
c
      integer local1, local2, u0
      integer sand, i, j
      sand(i,j)=iand(ishft(i,-j),1)
c
      do 20 i=1,((nuw+63)/64)
c
         u0 = (i-1)*64
         local1=p(2,i)
         do 11 j=31,0,-1
            u0=u0+1
            u(u0)=sand(local1,j)
11       continue
c
         u0 = (i-1)*64 + 32
         local2=p(1,i)
         do 12 j=31,0,-1
            u0=u0+1
            u(u0)=sand(local2,j)
12       continue
c
20    continue
#elif defined(FUJITSU_VP)
c
c  the following code may be used for 32-bit integer, big-endian
c   machines that use vax-type bit operators.
c  11-jun-89  written by ron shepard.
c
c  the general code breaks the alliant compiler. use this instead. -rls
      integer p(2,*)
      integer u(*)
c
      integer local1, local2, u0
      integer sand, i, j
      sand(i,j)=iand(ishft(i,-j),1)
c
      do 20 i=1,((nuw+63)/64)
c
         u0 = (i-1)*64
         local1=p(1,i)
         do 11 j=31,0,-1
            u0=u0+1
            u(u0)=sand(local1,j)
11       continue
c
         u0 = (i-1)*64 + 32
         local2=p(2,i)
         do 12 j=31,0,-1
            u0=u0+1
            u(u0)=sand(local2,j)
12       continue
c
20    continue
#else 
c  general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer u(*)
c
      integer i, ipower, j, u0
      integer local1, local2, local3, local4
      integer*2 llocl1(2), llocl2(2), llocl3(2), llocl4(2)
      equivalence (llocl1,local1),(llocl2,local2)
      equivalence (llocl3,local3),(llocl4,local4)
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
      intrinsic mod
c
      local1=0
      local2=0
      local3=0
      local4=0
c
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         u0=0
         do 10 i=1,((nuw+63)/64)
c
            llocl1(1)=p(4,i)
            ipower=2**15
            do 1 j=1,16
               u0=u0+1
               u(u0)=mod(local1/ipower,2)
               ipower=ipower/2
1           continue
c
            llocl2(1)=p(3,i)
            ipower=2**15
            do 2 j=1,16
               u0=u0+1
               u(u0)=mod(local2/ipower,2)
               ipower=ipower/2
2           continue
c
            llocl3(1)=p(2,i)
            ipower=2**15
            do 3 j=1,16
               u0=u0+1
               u(u0)=mod(local3/ipower,2)
               ipower=ipower/2
3           continue
c
            llocl4(1)=p(1,i)
            ipower=2**15
            do 4 j=1,16
               u0=u0+1
               u(u0)=mod(local4/ipower,2)
               ipower=ipower/2
4           continue
c
10       continue
      else
c        # ...big-endian.
         u0=0
         do 20 i=1,((nuw+63)/64)
c
            llocl1(2)=p(1,i)
            ipower=2**15
            do 11 j=1,16
               u0=u0+1
               u(u0)=mod(local1/ipower,2)
               ipower=ipower/2
11          continue
c
            llocl2(2)=p(2,i)
            ipower=2**15
            do 12 j=1,16
               u0=u0+1
               u(u0)=mod(local2/ipower,2)
               ipower=ipower/2
12           continue
c
            llocl3(2)=p(3,i)
            ipower=2**15
            do 13 j=1,16
               u0=u0+1
               u(u0)=mod(local3/ipower,2)
               ipower=ipower/2
13           continue
c
            llocl4(2)=p(4,i)
            ipower=2**15
            do 14 j=1,16
               u0=u0+1
               u(u0)=mod(local4/ipower,2)
               ipower=ipower/2
14           continue
c
20       continue
      endif
#endif 
      return
      end
