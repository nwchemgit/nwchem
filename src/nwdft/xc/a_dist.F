      Subroutine a_dist(xyz, Rij, nctrs,dopts)
c
C$Id: a_dist.F,v 1.12 2001-03-08 03:00:07 edo Exp $
c
      implicit none
c
#include "cdft.fh"
c
      integer nctrs
c
c     Cartesian Coordinates
c
      double precision xyz(3,nctrs)
c
c     Internuclear Distances
c
      double precision Rij((nctrs*(nctrs+1))/2)
c
      integer i, j, ij
      logical dopts
c
c     Evaluate the distance between each pair of centers.
c
c     pt charges are "tagged" with -1
c
      ij = 0
      do 20 i = 1, nctrs
         do 10 j = 1, i
            ij = ij + 1
            if (dopts.and.(iatype_pt_chg(i).or.iatype_pt_chg(j)))then
               Rij(ij) = -1.0d0
            elseif (i.eq.j)then
               Rij(ij) = 0.0d0
            else
               Rij(ij) = 1.0d0/sqrt(  (xyz(1,i) - xyz(1,j))**2 +
     &                                (xyz(2,i) - xyz(2,j))**2 +
     &                                (xyz(3,i) - xyz(3,j))**2   )
            endif
   10    continue
   20 continue
c
      return
      end
