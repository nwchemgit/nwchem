      Subroutine dft_gather(n,a,b,index,idiff)

C$Id: dft_gather.F,v 1.1 1995-04-27 03:16:42 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Dimension a(n),b(*),index(n)

c*******************************************************************************

      do 10 i = 1,n
       a(i) = b(index(i)-idiff)
   10 continue

      end
