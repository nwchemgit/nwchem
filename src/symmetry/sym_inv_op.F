      subroutine sym_inv_op(geom, op, opinv)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom, op, opinv
c
      opinv = sym_op_inv(op, geom)
c
      end
      subroutine sym_init_inv_op(geom)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom
      integer i
c
      do i = 1, sym_num_ops(geom)
         call sym_find_inv_op(geom, i, sym_op_inv(i,geom))
      enddo
c
      end
      subroutine sym_find_inv_op(geom, op, opinv)
C$Id$
      implicit none
#include "errquit.fh"
#include "sym.fh"
      integer geom, op, opinv
c
c     return in opinv the index of the operation that is
c     the inverse of op
c
      double precision r(3), rop(3), rinv(3), dist
      integer nops
c
      r(1) = 3.14159d0          ! A point not on a symmetry axis
      r(2) = 2.71828d0
      r(3) = 1.0d0
c
      nops = sym_number_ops(geom)
      call sym_apply_cart_op(geom, op, r, rop)
c
      do opinv = 1, nops
         call sym_apply_cart_op(geom, opinv, rop, rinv)
         dist = sqrt((r(1)-rinv(1))**2 + (r(2)-rinv(2))**2 + 
     $        (r(3)-rinv(3))**2)
         if (dist .lt. 1.0d-8) return
      enddo
c
c      call errquit('sym_inv_op: failed to find inverse for',op,
c     &       CALC_ERR)
c
      end
