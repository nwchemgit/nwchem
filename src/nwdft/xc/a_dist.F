      Subroutine a_dist(xyz, Rij, ncenters)

C$Id: a_dist.F,v 1.7 1995-11-03 00:07:59 og845 Exp $
      implicit none

#include "quad.fh"
      integer ncenters
c--> Cartesian Coordinates

      double precision xyz(3,ncenters)

c--> Internuclear Distances

      double precision Rij((ncenters*(ncenters+1))/2)

      integer i,j,ij
c
c Evaluate the distance between each pair of centers.
c
c***************************************************************************
c
c     eliminate pt charges from list
c
      
      i=0
  123 i=i+1

      if (aratio(iatype(i),1).lt.-1200.d0) then
        do j = i,ncenters-1
          xyz(1,j) = xyz(1,j+1)
          xyz(2,j) = xyz(2,j+1)
          xyz(3,j) = xyz(3,j+1)
        enddo
        ncenters=ncenters-1
        i=i-1
      endif
      if(i.lt.ncenters) goto 123
      
      ij = 0
      do 20 i = 1,ncenters
        do 10 j = 1,i
          ij = ij + 1
          if(i.eq.j)then
            Rij(ij) = 0.0d0
          else
            Rij(ij) = 1.0d0/sqrt(  (xyz(1,i) - xyz(1,j))**2 +
     &           (xyz(2,i) - xyz(2,j))**2 +
     &           (xyz(3,i) - xyz(3,j))**2   )
          endif
   10   continue
   20 continue
      return
      end
