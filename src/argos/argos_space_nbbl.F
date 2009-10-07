      subroutine argos_space_nbbl(n)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer n
c
      n=nbbl
      if(lqmd) n=0
c
      return
      end
