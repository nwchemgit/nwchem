      Subroutine a_dist(xyz, Rij, nctrs_no_pt_chgs)
c
C$Id: a_dist.F,v 1.8 1997-03-14 23:04:45 d3h449 Exp $
c
      implicit none
c
#include "cdft.fh"
c
      integer nctrs_no_pt_chgs
c
c---> Cartesian Coordinates
c
      double precision xyz(3,nctrs_no_pt_chgs)
c
c---> Internuclear Distances
c
      double precision Rij((nctrs_no_pt_chgs*(nctrs_no_pt_chgs+1))/2)
c
      integer i,j,ij
c
c     Evaluate the distance between each pair of centers.
c
c     eliminate pt charges from list
c
      i = 0
  123 i = i + 1
c
      if (aratio(iatype(i),1).lt.-1200.d0)then
         do j = i, nctrs_no_pt_chgs-1
            xyz(1,j) = xyz(1,j+1)
            xyz(2,j) = xyz(2,j+1)
            xyz(3,j) = xyz(3,j+1)
         enddo
         nctrs_no_pt_chgs = nctrs_no_pt_chgs - 1
         i = i - 1
      endif
c
      if (i.lt.nctrs_no_pt_chgs)goto 123
c      
      ij = 0
      do 20 i = 1, nctrs_no_pt_chgs
         do 10 j = 1, i
            ij = ij + 1
            if (i.eq.j)then
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
