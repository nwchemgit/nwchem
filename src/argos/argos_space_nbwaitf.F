      subroutine argos_space_nbwaitf()
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      call argos_space_nbwait3(int_mb(i_bb))
c
      return
      end
c $Id$
