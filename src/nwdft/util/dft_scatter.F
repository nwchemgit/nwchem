      Subroutine dft_scatter(n,a,index,b)

C$Id: dft_scatter.F,v 1.1 1995-04-27 03:16:43 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Dimension a(*),index(n),b(n)

c*******************************************************************************

      do 10 i = 1,n
       a(index(i)) = b(i)
   10 continue

      end
