      double precision function dabsmax(n, a)
C$Id: dabsmax.F,v 1.4 1997-06-03 21:45:37 d3e129 Exp $
      implicit none
c
      integer n
*     integer i 
      double precision a(n)
      integer idamax
      external idamax
c
*      dabsmax = 0.0d0
*      do i = 1, n
*         dabsmax = max(dabsmax, abs(a(i)))
*      enddo
c use idamax
c
      dabsmax = abs(a(idamax(n,a,1)))
c
      end
