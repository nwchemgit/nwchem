      subroutine plab32( p, u, nuw )
c
c  pack integral labels from u(*) into p(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+1)/2)*2 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  written by ron shepard.
c  version: 14-mar-89
c
      implicit none
      integer nuw
#if defined(CRAY)
      real*8 p(*)
      integer u(nuw)
c
      integer nuw32
c
      nuw32=((nuw+1)/2)*2
      if ( nuw32 .ne. 0 ) call pack(p,32,u,nuw32)

cevb Need to add definitions for oder, und, links and rechts
cevb for other platforms.  The bitwise functions are
cevb F90/F95 instrinsics

#elif defined(INT64)
      integer p(*), u(2,*)
      integer j,m32,i,nuw32
      parameter(m32=2**32-1)
      integer moder,oder,links,und
#  if defined(T3E64) || defined(T3D)
      moder(i,j)=or(i,and(j,m32))
      oder (i,j)=or(i,j)
      und  (i,j)=and(i,j)
      links (i,j)= shiftl(i,j)
#  elif defined(SGIPOWER) || defined(DECALPHA)
      moder(i,j)=ior(i,kiand(j,m32))
      oder (i,j)=ior(i,j)
      und  (i,j)=iand(i,j)
      links (i,j)=ishft(i,j)
#  elif defined(F90) || defined(F95)
      moder(i,j)=ior(i,iand(j,m32))
      oder (i,j)=ior(i,j)
      und  (i,j)=iand(i,j)
      links (i,j)=ishft(i,j)
#  endif 


      nuw32=(nuw+1)/2
      do 10 i=1,nuw32
      p(i)=moder(links(u(1,i),32), und(u(2,i),m32))
 10   continue
#else 
      integer p(2,*)
      integer u(2,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+1)/2)
            p(2,i)=u(1,i)
            p(1,i)=u(2,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+1)/2)
            p(1,i)=u(1,i)
            p(2,i)=u(2,i)
20       continue
      endif
#endif 
      return
      end
