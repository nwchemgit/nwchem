      subroutine argos_diana_rdssl()
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
c
      filsel=card(8:80)
      call argos_diana_readsel(lfnsel,filsel,int_mb(i_wsel),
     + int_mb(i_ssel))
      call argos_diana_wcount(int_mb(i_wsel),int_mb(i_ndxw))
c
      return
      end
c $Id$
