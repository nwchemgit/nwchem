      double precision function dabssum(n, a)
c $Id$
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
