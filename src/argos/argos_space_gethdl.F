      subroutine argos_space_gethdl(ibbl,lhandl,lself,
     + iwfr,iwto,isfr,isto)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      integer iwfr,iwto,isfr,isto
      logical lself
c
      call argos_space_gthdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto,
     + int_mb(i_bb),int_mb(i_ipl))
c
      return
      end
c $Id$
