      subroutine sym_apply_cart_op(geom, opnum, r, r1)
C$Id$
      implicit none
c
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom, opnum
      double precision r(3), r1(3), tmp(3), tmp2(3)
      integer i, j
c
c     Apply operator to vector of CARTESIAN coordinates
c
      do i = 1, 3
         tmp(i) = 0.0d0
         do j = 1, 3
            tmp(i) = tmp(i) + amatrix_inv(i,j,geom)*r(j)
         enddo
      enddo
      call sym_apply_op(geom, opnum, tmp, tmp2)
      do i = 1, 3
         r1(i) = 0.0d0
         do j = 1, 3
            r1(i) = r1(i) + amatrix(i,j,geom)*tmp2(j)
         enddo
      enddo
c
*      write(6,*) ' sym_apply_cart_op r ', r
*      write(6,*) ' sym_apply_cart_op op', opnum
*      write(6,*) ' sym_apply_cart_op r1', r1
c
      end
