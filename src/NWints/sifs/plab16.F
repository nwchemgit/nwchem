      subroutine plab16( p, u, nuw )
c
c  pack integral labels from u(*) into p(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+3)/4)*4 ) are referenced.
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
      integer nuw16
c
      nuw16=((nuw+3)/4)*4
      if ( nuw16 .ne. 0 ) call pack(p,16,u,nuw16)

cevb Need to add definitions for oder, und, links and rechts
cevb for other platforms.  The bitwise functions are
cevb F90/F95 instrinsics

#elif defined(INT64)
      integer p(*), u(4,*)
      integer m16,m32
      parameter (m16=2**16-1, m32=2**32-1)
      integer i,j,jscr32a, jscr32b,nuw16
      integer moder,oder,links,und
#  if defined(T3E64) || defined(T3D)
      moder(i,j)=or(i,and(j,m16))
      oder (i,j)=or(i,j)
      und  (i,j)=and(i,j)
      links (i,j)= shiftl(i,j)
#  elif defined(SGIPOWER) || defined(DECALPHA)
      moder(i,j)=ior(i,kiand(j,m16))
      oder (i,j)=ior(i,j)
      und  (i,j)=iand(i,j)
      links (i,j)=ishft(i,j)
#  elif defined(F90) || defined(F95)
      moder(i,j)=ior(i,iand(j,m16))
      oder (i,j)=ior(i,j)
      und  (i,j)=iand(i,j)
      links (i,j)=ishft(i,j)
#  endif 

       nuw16=(nuw+3)/4
        do i=5 - nuw16*4 +nuw,4
        u(i,nuw16)=0
         enddo
       do 10 i=1,nuw16
       jscr32a= moder(links(u(1,i),16),u(2,i))
       jscr32b= moder(links(u(3,i),16),u(4,i))
       p(i)= oder(links(jscr32a,32),und(jscr32b,m32) )
  10   continue
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
      longw = 1
      if ( shortw(1) .eq. 1 ) then
c        # ...little-endian.
         do 10 i=1,((nuw+3)/4)
            p(4,i)=u(1,1,i)
            p(3,i)=u(1,2,i)
            p(2,i)=u(1,3,i)
            p(1,i)=u(1,4,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+3)/4)
            p(1,i)=u(2,1,i)
            p(2,i)=u(2,2,i)
            p(3,i)=u(2,3,i)
            p(4,i)=u(2,4,i)
20       continue
      endif
#endif 
      return
      end
