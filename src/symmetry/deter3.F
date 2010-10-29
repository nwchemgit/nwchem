      double precision function deter3(r)
C$Id$
      implicit none
      double precision r(3,3)
c
c     return the determinant of a 3x3 matrix
c
      deter3 = 
     $     r(1,1)*(r(2,2)*r(3,3)-r(2,3)*r(3,2)) -
     $     r(1,2)*(r(2,1)*r(3,3)-r(2,3)*r(3,1)) + 
     $     r(1,3)*(r(2,1)*r(3,2)-r(2,2)*r(3,1))
c
      end
