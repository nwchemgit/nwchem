      Subroutine a_dist(xyz, Rij, ncenter)

C$Id: a_dist.F,v 1.5 1995-04-26 18:12:38 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)


c--> Cartesian Coordinates

      Dimension xyz(3,ncenter)

c--> Internuclear Distances

      Dimension Rij((ncenter*(ncenter+1))/2)
c
c Evaluate the distance between each pair of centers.
c
c***************************************************************************
      ij = 0
      do 20 i = 1,ncenter
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
