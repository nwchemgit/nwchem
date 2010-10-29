      Subroutine dft_scatter(n,a,index,b)

C$Id$
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
