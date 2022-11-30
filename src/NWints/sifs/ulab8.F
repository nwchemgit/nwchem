      subroutine ulab8( p, u, nuw )
c
c  unpack integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u(1: ((nuw+7)/8)*8 ) are referenced.
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
      integer nuw8
c
      nuw8=((nuw+7)/8)*8
      if ( nuw8 .ne. 0 ) call unpack(p,8,u,nuw8)

cevb Need to add definitions for oder, und, links and rechts
cevb for other platforms.  The bitwise functions are
cevb F90/F95 instrinsics

#elif defined(INT64)
      integer p(*)
      integer u(8,*)
c
      integer nuw8
c
      integer   m8
      parameter(m8=2**8-1)
c
       integer und,rechts
        integer i,j
#  if defined(T3E64) || defined(T3D)
       und(i,j) = and(i,j)
       rechts(i,j) = shiftr(i,j)
#  elif defined(SGIPOWER) || defined(DECALPHA)
       und(i,j) = iand(i,j)
       rechts(i,j)=rshift(i,j)
#  elif defined(F90) || defined(F95)
       und(i,j) = iand(i,j)
       rechts(i,j)=ishft(i,-j)
#  endif 
cgk debug
*     write(*,*)'gk: in ulab8 int64 section'
cgk end

c
c     call izero(nuw,u,1)
      nuw8=(nuw+7)/8
      do 10 i=1,nuw8
         u(1,i)=und(rechts(p(i),56),m8)
         u(2,i)=und(rechts(p(i),48),m8)
         u(3,i)=und(rechts(p(i),40),m8)
         u(4,i)=und(rechts(p(i),32),m8)
         u(5,i)=und(rechts(p(i),24),m8)
         u(6,i)=und(rechts(p(i),16),m8)
         u(7,i)=und(rechts(p(i),8),m8)
         u(8,i)=und(       p(i),    m8)
10    continue
#elif defined(OLDOLDSUN)
      integer p(2,*)
      integer u(8,*)
c
      integer i, nuw8
c
      integer   m8
      parameter(m8=2**8-1)
c
      integer   and, rshift
      intrinsic and, rshift
c
      nuw8=(nuw+7)/8
      do 10 i=1,nuw8
         u(1,i)=and(rshift(p(1,i),24),m8)
         u(2,i)=and(rshift(p(1,i),16),m8)
         u(3,i)=and(rshift(p(1,i), 8),m8)
         u(4,i)=and(       p(1,i),    m8)
         u(5,i)=and(rshift(p(2,i),24),m8)
         u(6,i)=and(rshift(p(2,i),16),m8)
         u(7,i)=and(rshift(p(2,i), 8),m8)
         u(8,i)=and(       p(2,i),    m8)
10    continue
#else 
c  general byte-addressable 32-bit integer machines.
      logical*1 p(8,*)
      logical*1 u(4,8,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
cgk debug
*     write(*,*)'gk: in ulab8 general byte-addressable section'
cgk end
      call izero(nuw,u,1)
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+7)/8)
            u(1,1,i)=p(8,i)
            u(1,2,i)=p(7,i)
            u(1,3,i)=p(6,i)
            u(1,4,i)=p(5,i)
            u(1,5,i)=p(4,i)
            u(1,6,i)=p(3,i)
            u(1,7,i)=p(2,i)
            u(1,8,i)=p(1,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+7)/8)
            u(4,1,i)=p(1,i)
            u(4,2,i)=p(2,i)
            u(4,3,i)=p(3,i)
            u(4,4,i)=p(4,i)
            u(4,5,i)=p(5,i)
            u(4,6,i)=p(6,i)
            u(4,7,i)=p(7,i)
            u(4,8,i)=p(8,i)
20       continue
      endif
#endif 
      return
      end
