      subroutine ulab32( p, u, nuw )
c
c  unpack integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+1)/2)*2 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  20-jul-91 cray nuw check added. -galen gawboy/rls
c  written by ron shepard.
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
      if ( nuw32 .ne. 0 ) call unpack(p,32,u,nuw32)

cevb Need to add definitions for und and rechts
cevb for other platforms.  The bitwise functions are
cevb F90/F95 instrinsics

#elif defined(INT64)
       integer  p(*), u(2,*)
       integer m32,nuw32
       parameter (m32=2**32-1)
       integer und,rechts
        integer i,j
#  if defined(T3E64) || defined(T3D)
       und(i,j) = and(i,j)
       rechts(i,j) = shiftr(i,j)
#  elif defined(SGIPOWER) || defined(DECALPHA)
       und(i,j) = iand(i,j)
       rechts(i,j)=ishft(i,-j)
#  elif defined(F90) || defined(F95)
       und(i,j) = iand(i,j)
       rechts(i,j)=ishft(i,-j)
#  endif 
       nuw32=((nuw+1)/2)
       call izero(nuw,u,1)
       do 10 i=1, nuw32
         u(1,i) = und(rechts(p(i),32),m32)
         u(2,i) = und(p(i),m32)
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
            u(2,i)=p(1,i)
            u(1,i)=p(2,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+1)/2)
            u(1,i)=p(1,i)
            u(2,i)=p(2,i)
20       continue
      endif
#endif 
      return
      end
