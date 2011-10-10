      subroutine argos_diana_setiomode()
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      read(card(8:12),1000) iomode
 1000 format(i5)
c
      call ga_brdcst(mag_d01,iomode,ma_sizeof(mt_int,1,mt_byte),0)
c
      return
      end
c $Id$
