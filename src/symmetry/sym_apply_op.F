      subroutine sym_apply_op(geom, opnum, r, r1)
C$Id: sym_apply_op.F,v 1.2 1995-02-02 23:23:18 d3g681 Exp $
      implicit none
c
#include "mafdecls.h"
#include "geomP.fh"
c
      integer geom, opnum
      double precision r(3), r1(3)
      integer i, j
c
c     Apply operator to vector of FRACTIONAL coordinates
c
      do i = 1, 3
        r1(i) = sym_ops((opnum-1)*3+i,4,geom)
        do j = 1, 3
          r1(i) = r1(i) + sym_ops((opnum-1)*3+i,j,geom)*r(j)
        enddo
      enddo
c
      write(6,*) ' sym_apply_op r ', r
      write(6,*) ' sym_apply_op op', opnum
      write(6,*) ' sym_apply_op r1', r1
c
      end
