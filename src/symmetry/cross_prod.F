      subroutine cross_product(a,b,c)
C$Id$
      implicit none
      double precision a(3), b(3) ! [input]
      double precision c(3)     ! [output]
c
c     return c = a ^ b ... i.e. the cross product of 
c     a and b.
c
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = -(a(1)*b(3) - a(3)*b(1))
      c(3) = a(1)*b(2) - a(2)*b(1)
c
      end
