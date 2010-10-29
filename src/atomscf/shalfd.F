      subroutine shalfd(s,v,n)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c.......................................................................
c
c     make s**(-1/2) for transformation to orthonormal basis.
c.......................................................................
      integer n
      real *8 s(*), v(n,n)
c
c      real *8  scr1(50), scr2(50)
      integer k_scr1,l_scr1,k_scr2,l_scr2
      real *8 hlp
      integer icount, nn, nnd
      integer i, j, k
c
      if(.not.ma_push_get(MT_DBL,n,'scr1',l_scr1,k_scr1))
     A call errquit(' shalfd: mapushget failed ',0,0)
      if(.not.ma_push_get(MT_DBL,n,'scr1',l_scr2,k_scr2))
     A call errquit(' shalfd: mapushget failed ',0,0)
      icount = 0
      nnd = n*(n+1)/2
      nn = n*n
c      if (n.gt.50)
c     &    call errquit('atomscf:dimensioning error in shalfd',n,
c     &       UNKNOWN_ERR)
      do i = 1 , n
         do j = 1 , i
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
         enddo
         v(i,i) = 1.0d0
      enddo
      call jacod(s,v,n,nnd,nn,1,n,
     D     dbl_mb(k_scr1),dbl_mb(k_scr2),n)
      icount = 1
      do i = 1 , n
         dbl_mb(k_scr1+i-1) = 0d0
         if(s(icount).gt.0d0)
     D   dbl_mb(k_scr1+i-1) = 1.0d0/sqrt(s(icount))
         icount = icount + i + 1
      enddo
      icount = 0
      do i = 1 , n
         do j = 1 , i
            icount = icount + 1
            hlp = 0.0d0
            do k = 1 , n
               hlp = hlp + v(i,k)*dbl_mb(k_scr1+k-1)*v(j,k)
            enddo
            s(icount) = hlp
         enddo
      enddo
      icount = 0
      do i = 1 , n
         do j = 1 , i
            icount = icount + 1
            v(j,i) = s(icount)
            v(i,j) = s(icount)
         enddo
      enddo
      if(.not.ma_chop_stack(l_scr1))
     A call errquit(' shalfd: machopstack failed ',0,0)
      return
      end
