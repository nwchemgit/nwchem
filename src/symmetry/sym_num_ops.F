      integer function sym_number_ops(geom)
C$Id: sym_num_ops.F,v 1.2 1999-07-14 00:29:11 d3e129 Exp $
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom
c
      sym_number_ops = sym_num_ops(geom)
c
      end
