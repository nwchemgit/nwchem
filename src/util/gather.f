
      subroutine gather(n,a,b,indx)
*
* $Id: gather.f,v 1.2 1997-10-31 20:45:32 d3e129 Exp $
*
      integer n, indx(n)
      double precision a(n),b(*)
      integer i
      
      do i=1,n
        a(i) = b(indx(i))
      enddo
      return
      end

      






