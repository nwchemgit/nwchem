      double precision function dabsmax(n, a)
C$Id: dabsmax.F,v 1.2 1997-05-29 00:49:09 d3e129 Exp $
      implicit none
c
      integer n, i
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
      dabsmax = a(idamax(n,a,1))
c
      end
