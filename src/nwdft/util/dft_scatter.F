      Subroutine dft_scatter(n,a,index,b)

C$Id: dft_scatter.F,v 1.2 1995-12-22 22:45:16 og845 Exp $
      implicit none
      
      integer n
      integer index(n)
      double precision a(*),b(n)
      integer i

c*******************************************************************************

      do 10 i = 1,n
       a(index(i)) = b(i)
   10 continue

      end
