      double precision function dabssum(n, a)
c $Id: dabssum.F,v 1.1 1995-04-11 19:00:10 rg240 Exp $
      implicit none
c
      integer n, i
      double precision a(n)
c
      dabssum = 0.0d0
      do i = 1, n
         dabssum = dabssum + abs(a(i))
      enddo
c
      end
