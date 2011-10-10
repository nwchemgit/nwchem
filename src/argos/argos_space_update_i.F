      subroutine argos_space_update_i(numsa,isl,numwm,iwl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer numsa,isl(msa,mis2)
      integer numwm,iwl(mwm,miw2)
c
      call argos_space_upd_i(numsa,isl,int_mb(i_pack),numwm,
     + iwl,int_mb(i_packw))
c
      return
      end
c $Id$
