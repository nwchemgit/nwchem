      Subroutine a_dist(xyz, Rij, nctrs,dopts)
c
C$Id$
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
      rij(1)=0d0
      if(nctrs.eq.1) return
      if (dopts.and.iatype_pt_chg(1)) rij(1) = -1.0d0
      ij = 1
      if(dopts) then
         do i = 2, nctrs
            do  j = 1, i-1
                ij = ij + 1
               if (iatype_pt_chg(i).or.iatype_pt_chg(j))then
                  Rij(ij) = -1.0d0
               else
                  Rij(ij) = 1.0d0/sqrt(  (xyz(1,i) - xyz(1,j))**2 +
     &                 (xyz(2,i) - xyz(2,j))**2 +
     &                                (xyz(3,i) - xyz(3,j))**2   )
               endif
            enddo
            ij=ij+1
            Rij(ij) = 0.0d0
         enddo
      else
         do i = 2, nctrs
!DEC$ LOOP COUNT(100)
            do  j = 1, i-1
               ij = ij + 1
               Rij(ij) = 1.0d0/sqrt(
     &              (xyz(1,i) - xyz(1,j))*(xyz(1,i) - xyz(1,j)) +
     &              (xyz(2,i) - xyz(2,j))*(xyz(2,i) - xyz(2,j)) +
     &              (xyz(3,i) - xyz(3,j))*(xyz(3,i) - xyz(3,j)))
            enddo
            ij=ij+1
            Rij(ij) = 0.0d0
         enddo
      endif
c
      return
      end
