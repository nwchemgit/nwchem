      Subroutine a_dist(xyz, Rij, ncenter)

C$Id: a_dist.F,v 1.3 1995-02-02 23:20:22 d3g681 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)


c--> Cartesian Coordinates

      Dimension xyz(3,ncenter)

c--> Internuclear Distances

      Dimension Rij((ncenter*(ncenter-1))/2)
      Logical LResult
c
c Evaluate the distance between each pair of centers.
c
c***************************************************************************
      ij = 0
      do 10 j = 2,ncenter
      do 10 i = 1,j-1
       ij = ij + 1

        Rij(ij) = sqrt(  (xyz(1,i) - xyz(1,j))**2 +
     &                   (xyz(2,i) - xyz(2,j))**2 +
     &                   (xyz(3,i) - xyz(3,j))**2   )

   10 continue

      end
