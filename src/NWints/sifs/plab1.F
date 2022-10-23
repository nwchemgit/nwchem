      subroutine plab1( p, u, nuw )
c
c  pack a bit vector.
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
      if ( nuw64 .ne. 0 ) call pack(p,1,u,nuw64)

cevb Need to add definitions for oder, und, links and rechts
cevb for other platforms.  The bitwise functions are
cevb F90/F95 instrinsics

#elif defined(INT64)
        integer p(*), u(64,*)
        integer oder,local,i,j,mpack,nuw64
#  if defined(T3E64) || defined(T3D)
        mpack(i,j)=shiftl(and(i,1),j)
        oder(i,j) = or(i,j)
#  elif defined(SGIPOWER) || defined(DECALPHA)
        mpack(i,j) =ishft(iand(i,1),j)
        oder(i,j)  = ior(i,j)
#  elif defined(F90) || defined(F95)
        mpack(i,j) = ishft(iand(i,1),j)
        oder(i,j)  = ior(i,j)
#  endif 
        nuw64= (nuw+63)/64
        do 10 i=1,nuw64
         p(i)=0
        do 11 j=63,0,-1
         p(i)= or(p(i), mpack(u(64-j,i),j))
   11   continue
   10   continue
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
      integer sor, i, j
      sor(i,j) = ior(ishft(i,1), iand(j,1) )
c
      do 20 i=1,((nuw+63)/64)
c
         u0 = (i-1)*64
         local1 = 0
         do 11 j = 1, 32
            local1 = sor(local1,u(u0+j))
11       continue
         p(2,i) = local1
c
         u0 = (i-1)*64 + 32
         local2 = 0
         do 12 j = 1, 32
            local2 = sor(local2,u(u0+j))
12       continue
         p(1,i) = local2
c
20    continue
#elif defined(FUJITSU_VP)
c
c  the following code may be used for 32-bit integer, big-endian
c  machines that use vax-type bit operators.
c  30-nov-89  written by ron shepard.
c
c  the general code breaks the titan compiler. use this version
c  instead. -rls
      integer p(2,*)
      integer u(*)
c
      integer local1, local2, u0
      integer sor, i, j
      sor(i,j) = ior(ishft(i,1), iand(j,1) )
c
      do 20 i=1,((nuw+63)/64)
c
         u0 = (i-1)*64
         local1 = 0
         do 11 j = 1, 32
            local1 = sor(local1,u(u0+j))
11       continue
         p(1,i) = local1
c
         u0 = (i-1)*64 + 32
         local2 = 0
         do 12 j = 1, 32
            local2 = sor(local2,u(u0+j))
12       continue
         p(2,i) = local2
c
20    continue
#else 
c  general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer u(*)
c
      integer i, u0
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
c  define bit-packing functions recursively for maximum pipelining...
      intrinsic mod
      integer lp2,lp4,lp8,lp16
      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16
      lp2(i1,i2) = mod(i1,2) * 2 + mod(i2,2)
      lp4(i1,i2,i3,i4) = (lp2(i1,i2) * 4 + lp2(i3,i4))
      lp8(i1,i2,i3,i4,i5,i6,i7,i8)=
     & (lp4(i1,i2,i3,i4) * 16 + lp4(i5,i6,i7,i8))
      lp16(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16)=
     & ( lp8(i1,i2,i3,i4,i5,i6,i7,i8) * 256 +
     & lp8(i9,i10,i11,i12,i13,i14,i15,i16) )
c
      longw = 1
      if(shortw(1).eq.1)then
c        # ...little-endian.
         do 10 i=1,((nuw+63)/64)
c
            u0=(i-1)*64
            local1=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(4,i)=llocl1(1)
c
            u0=(i-1)*64+16
            local2=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(3,i)=llocl2(1)
c
            u0=(i-1)*64+32
            local3=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(2,i)=llocl3(1)
c
            u0=(i-1)*64+48
            local4=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(1,i)=llocl4(1)
c
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+63)/64)
c
            u0=(i-1)*64
            local1=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(1,i)=llocl1(2)
c
            u0=(i-1)*64+16
            local2=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(2,i)=llocl2(2)
c
            u0=(i-1)*64+32
            local3=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(3,i)=llocl3(2)
c
            u0=(i-1)*64+48
            local4=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
            p(4,i)=llocl4(2)
c
20       continue
      endif
#endif 
      return
      end
