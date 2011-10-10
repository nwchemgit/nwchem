      subroutine argos_diana_savsel()
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
c
      filsel=card(8:80)
      call argos_diana_wrtsel(lfnsel,filsel,int_mb(i_wsel),
     + int_mb(i_isel))
c
      return
      end
c $Id$
