      Subroutine grid_neighb_list(x, xyzm,  nctrs_scr,
     ,     ictr,indx)
c
C$Id$
c
c     xc_neighb_list removes point charges and centers far away
c     from the current integration center.
c
      implicit none
c
#include "cdft.fh"
c
      integer nctrs_scr
      integer indx(*) ! [in/out]
c
c     Cartesian Coordinates
c
      double precision x(3)
      double precision xyzm(3, *)
      integer ictr
      integer i, i_a, j
      double precision dist
      logical ok
c
c     eliminate distant centers and point charges
c
      i = 0
      i_a = 0
  123 i = i + 1
      i_a = i_a + 1
      ok = iatype_pt_chg(i_a)
      dist = 1.d6
      if (.not.ok)
     &   dist = sqrt( (xyzm(1,i) - x(1))**2 +
     &                (xyzm(2,i) - x(2))**2 +
     &                (xyzm(3,i) - x(3))**2   )
      ok = ok.or.(dist.gt.radius)
c
c     get center 
c
      if (ok)then
         do j = i, nctrs_scr - 1
            xyzm(1,j) = xyzm(1,j+1)
            xyzm(2,j) = xyzm(2,j+1)
            xyzm(3,j) = xyzm(3,j+1)
            indx(j)=indx(j+1)
         enddo
         nctrs_scr = nctrs_scr - 1
         i = i - 1
      elseif (dist.lt.1d-3) then
         ictr=i
      endif
      if(i.lt.nctrs_scr) goto 123
      return
      end
