      double precision function dabsmax(n, a)
C$Id$
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
