      subroutine sym_ops_get(geom,numops,symops)
C$Id: sym_ops_get.F,v 1.2 1995-02-02 23:23:26 d3g681 Exp $
      implicit none
c
#include "mafdecls.h"
#include "geomP.fh"
c
      integer i,j,geom,numops,l,k
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
