      subroutine sym_ops_get(geom,numops,symops)
C$Id: sym_ops_get.F,v 1.3 1995-10-10 06:17:46 d3g681 Exp $
      implicit none
c
#include "mafdecls.h"
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
