      subroutine scatter(n,a,indx,b)
*
* $Id: scatter.f,v 1.2 1997-10-31 20:45:34 d3e129 Exp $
*
      integer n, indx(n)
      double precision a(*), b(n)
      integer i
      
      do i=1,n
        a(indx(i)) = b(i)
      enddo
      return
      end
