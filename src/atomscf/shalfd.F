      subroutine shalfd(s,v,n)
C$Id: shalfd.F,v 1.2 1995-02-02 23:09:10 d3g681 Exp $
      implicit none
c.......................................................................
c
c     make s**(-1/2) for transformation to orthonormal basis.
c.......................................................................
      integer n
      real *8 s(*), v(n,n)
c
      real *8  scr1(50), scr2(50)
      real *8 hlp
      integer icount, nn, nnd
      integer i, j, k
c
      icount = 0
      nnd = n*(n+1)/2
      nn = n*n
      if (n.gt.50) call errquit('dimensioning error in shalfd',0)
      do i = 1 , n
         do j = 1 , i
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
         enddo
         v(i,i) = 1.0d0
      enddo
      call jacod(s,v,n,nnd,nn,1,n,scr1,scr2,n)
      icount = 1
      do i = 1 , n
         scr1(i) = 1.0d0/sqrt(s(icount))
         icount = icount + i + 1
      enddo
      icount = 0
      do i = 1 , n
         do j = 1 , i
            icount = icount + 1
            hlp = 0.0d0
            do k = 1 , n
               hlp = hlp + v(i,k)*scr1(k)*v(j,k)
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
      return
      end
