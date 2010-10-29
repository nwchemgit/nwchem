      Subroutine dft_gather(n,a,b,index,idiff)

C$Id$
      implicit none
      
      integer n
      double precision a(n),b(*)
      integer index(n),idiff
      integer i

c****************************************************************************

      do 10 i = 1,n
       a(i) = b(index(i)-idiff)
   10 continue

      end
