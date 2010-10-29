      integer function sym_number_ops(geom)
C$Id$
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
      integer geom
c
      sym_number_ops = sym_num_ops(geom)
c
      end
