      subroutine ulab16( p, u, nuw )
c
c  unpack one-electron integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u(1: ((nuw+3)/4)*4 ) are referenced.
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
      integer nuw16
c
      nuw16=((nuw+3)/4)*4
      if ( nuw16 .ne. 0 ) call unpack(p,16,u,nuw16)

cevb Need to add definitions for oder, und, links and rechts
cevb for other platforms.  The bitwise functions are
cevb F90/F95 instrinsics

#elif defined(INT64)
       integer p(*), u(4,*)
       integer m16,nuw16
       parameter (m16 = 2**16-1)
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

       nuw16=((nuw+3)/4)
c      call izero (nuw , u, 1)
       do 10 i=1, nuw16
        u(1,i)=und(rechts(p(i),48),m16)
        u(2,i)=und(rechts(p(i),32),m16)
        u(3,i)=und(rechts(p(i),16),m16)
        u(4,i)=und(p(i),m16)
 10    continue
#else 
c     # general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer*2 u(2,4,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer longw
      integer*2 shortw(2)
      equivalence (longw,shortw)
c
      call izero ( nuw, u, 1 )
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+3)/4)
            u(1,1,i)=p(4,i)
            u(1,2,i)=p(3,i)
            u(1,3,i)=p(2,i)
            u(1,4,i)=p(1,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+3)/4)
            u(2,1,i)=p(1,i)
            u(2,2,i)=p(2,i)
            u(2,3,i)=p(3,i)
            u(2,4,i)=p(4,i)
20       continue
      endif
#endif 
      return
      end
