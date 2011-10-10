      subroutine argos_cafe_print_distar(lfnout)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer lfnout
c
      call argos_cafe_print_dstr(lfnout,byte_mb(i_snam),
     + int_mb(i_idrs),dbl_mb(i_rdrs))
c
      return
      end
c $Id$
