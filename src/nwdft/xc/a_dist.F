      Subroutine a_dist(geom, xyz, Rij, ncenter, charge, tags )

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

#include <geom.fh>

c****
c**** nwchem geometry handle
c****
      integer geom

c--> Cartesian Coordinates

      Dimension xyz(3,ncenter), charge(ncenter)

      character*16 tags(ncenter)

c--> Internuclear Distances

      Dimension Rij((ncenter*(ncenter-1))/2)
      Logical LResult
c
c Evaluate the distance between each pair of centers.
c
c***************************************************************************
      LResult = geom_cart_get( geom, ncenter, tags, xyz, charge )

      ij = 0
      do 10 j = 2,ncenter
      do 10 i = 1,j-1
       ij = ij + 1

        Rij(ij) = sqrt(  (xyz(1,i) - xyz(1,j))**2 +
     &                   (xyz(2,i) - xyz(2,j))**2 +
     &                   (xyz(3,i) - xyz(3,j))**2   )

   10 continue

      end
