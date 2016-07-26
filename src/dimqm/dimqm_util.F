c     
c     $Id: dft_scf.F 28382 2016-07-21 01:23:33Z niri $
c     
      subroutine dimqm_used(ldimqm_out)
      implicit none
      logical ldimqm_out
#include "dimqm.fh"
      ldimqm_out=ldimqm
      return
      end
      subroutine dimqm_getgo(dimqm_go_out)
      implicit none
      logical dimqm_go_out
#include "dimqm.fh"
      dimqm_go_out=dimqm_go
      return
      end
      subroutine dimqm_setgo(dimqm_go_in)
      implicit none
      logical dimqm_go_in
#include "dimqm.fh"
      dimqm_go=dimqm_go_in
      return
      end
      subroutine dimqm_energy(edimqm_out)
      implicit none
      double precision edimqm_out
#include "dimqm.fh"
      edimqm_out=edimqm
      return
      end
