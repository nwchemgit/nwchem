
      subroutine gather(n,a,b,indx)
*
* $Id$
*
      integer n, indx(n)
      double precision a(n),b(*)
      integer i
      
      do i=1,n
        a(i) = b(indx(i))
      enddo
      return
      end

      






