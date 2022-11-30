      subroutine plab8( p, u, nuw )
c
c  pack integral labels from u(*) into p(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+7)/8)*8 ) are referenced.
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
      integer nuw8
c
      nuw8=((nuw+7)/8)*8
      if ( nuw8 .ne. 0 ) call pack(p,8,u,nuw8)

cevb Need to add definitions for oder, und, links and rechts
cevb for other platforms.  The bitwise functions are
cevb F90/F95 instrinsics

#elif defined(INT64)
      integer p(*),m8,m32
      integer u(8,*)
      integer nuw8
c
      parameter( m32=2**32-1)
      parameter( m8=2**8-1)
c
C
      integer jscr8,jscr16,jscr16b
      integer moder,oder,links,und
        integer i,j
#  if defined(T3E64) || defined(T3D)
      moder(i,j)=or(i,and(j,m8))
      oder (i,j)=or(i,j)
      und  (i,j)=and(i,j)
      links (i,j)= shiftl(i,j)
#  elif defined(SGIPOWER) || defined(DECALPHA)
      moder(i,j)=ior(i,kiand(j,m8))
      oder (i,j)=ior(i,j)
      und  (i,j)=iand(i,j)
      links (i,j)=ishft(i,j)
#  elif defined(F90) || defined(F95)
      moder(i,j)=ior(i,iand(j,m8))
      oder (i,j)=ior(i,j)
      und  (i,j)=iand(i,j)
      links (i,j)=ishft(i,j)
#  endif 
c
      nuw8=(nuw+7)/8
       do  i=9-nuw8*8+nuw, 8
       u(i,nuw8)=0
       enddo
      do 10 i=1,nuw8
       jscr8=moder( links( u(1,i),8),u(2,i))
       jscr16=oder (links(jscr8,16), oder(links(u(3,i),8),u(4,i) ))
      jscr8= moder(links( u(5,i),8),u(6,i))
       jscr16b= oder(links(jscr8,16), moder(links(u(7,i),8),u(8,i)) )
       p(i) = oder(links(jscr16,32), und(jscr16b,m32) )

10    continue
#elif defined(OLDOLDSUN)
c  17-apr-91 mask added to fix garbage overwrite bug. -rls/tom kovar
      integer p(2,*)
      integer u(8,*)
c
      integer nuw8
c
      integer   or, lshift, and
      intrinsic or, lshift, and
c
      integer    m8
      parameter( m8=2**8-1 )
c
      integer i, j
      integer mor
      mor(i,j) = or( i, and(j,m8) )
c
      nuw8=(nuw+7)/8
      do 10 i=1,nuw8
         p(1,i)=mor(lshift(mor(lshift(mor(lshift(
     &    u(1,i),8),u(2,i)),8),u(3,i)),8),u(4,i))
         p(2,i)=mor(lshift(mor(lshift(mor(lshift(
     &    u(5,i),8),u(6,i)),8),u(7,i)),8),u(8,i))
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
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+7)/8)
            p(8,i)=u(1,1,i)
            p(7,i)=u(1,2,i)
            p(6,i)=u(1,3,i)
            p(5,i)=u(1,4,i)
            p(4,i)=u(1,5,i)
            p(3,i)=u(1,6,i)
            p(2,i)=u(1,7,i)
            p(1,i)=u(1,8,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+7)/8)
            p(1,i)=u(4,1,i)
            p(2,i)=u(4,2,i)
            p(3,i)=u(4,3,i)
            p(4,i)=u(4,4,i)
            p(5,i)=u(4,5,i)
            p(6,i)=u(4,6,i)
            p(7,i)=u(4,7,i)
            p(8,i)=u(4,8,i)
20       continue
      endif
#endif 
      return
      end
