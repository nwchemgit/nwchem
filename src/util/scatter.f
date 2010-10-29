      subroutine scatter(n,a,indx,b)
*
* $Id$
*
      integer n, indx(n)
      double precision a(*), b(n)
      integer i
      
      do i=1,n
        a(indx(i)) = b(i)
      enddo
      return
      end
