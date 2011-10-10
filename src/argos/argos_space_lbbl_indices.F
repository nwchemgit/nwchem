      subroutine argos_space_lbbl_indices()
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      call argos_space_lbbl_ind(int_mb(i_bb),int_mb(i_ipl),
     + int_mb(i_jpl))
c
      return
      end
c $Id$
