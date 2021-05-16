      subroutine sym_op_get(geom,opnum,matrix)
C$Id$
      implicit none
c
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer i,j,geom,opnum,opoff
c
      double precision matrix(3,4)
c
      opoff = (opnum-1)*3
      do j = 1, 4
         do i = 1, 3
            matrix(i,j) = sym_ops(opoff+i,j,geom)
         enddo
      enddo

      return
      end
