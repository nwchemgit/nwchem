      double precision function dabsmax(n, a)
C$Id: dabsmax.F,v 1.1 1995-04-11 19:00:09 rg240 Exp $
      implicit none
c
      integer n, i
      double precision a(n)
c
      dabsmax = 0.0d0
      do i = 1, n
         dabsmax = max(dabsmax, abs(a(i)))
      enddo
c
      end
