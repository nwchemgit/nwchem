      subroutine sym_ops_get(geom,numops,symops)
C$Id: sym_ops_get.F,v 1.5 1999-07-14 00:29:11 d3e129 Exp $
      implicit none
c
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer i,j,geom,numops
c
      double precision symops(numops*3,4)
c
      do i=1,numops*3
        do j=1,4
          symops(i,j)=sym_ops(i,j,geom)
        enddo
      enddo
c
      return
      end
