      subroutine tramad(a,b,c,mdima,mdim,scr)
C$Id$
      implicit none
c.......................................................................
c
c     this routine tranforms a to (b+)ab, where a is a lower
c     triangular symmetric matrix, and b is a square matrix.
c     the transformed matrix is returned in a.
c.......................................................................
      integer mdima, mdim
      real *8 a(mdima), b(mdim,mdim), c(mdima), scr(mdim)
c
      integer i, j, k, l, maxjl
      real *8 sum
c
      do i = 1 , mdim
c.......................................................................
c
c     generate i'th column of ab.
c.......................................................................
         do j = 1 , mdim
            sum = 0.0d0
            do l = 1 , mdim
               maxjl = max(j,l)
               sum = sum + a(maxjl*(maxjl-3)/2 + j + l)*b(l,i)
            enddo
            scr(j) = sum
         enddo
c.......................................................................
c
c     multiply this by rows of b+
c.......................................................................
         do j = 1 , i
            sum = 0.0d0
            do k = 1 , mdim
               sum = sum + scr(k)*b(k,j)
            enddo
            c(i*(i-1)/2+j) = sum
         enddo
      enddo
c.......................................................................
c
c     transfer to a for return.
c.......................................................................
      do i = 1 , mdima
         a(i) = c(i)
      enddo
      return
      end
