      subroutine argos_diana_tag(sgmnam,x)
c
c $Id: argos_diana_tag.F,v 1.2 2007-01-27 02:22:35 d3j191 Exp $
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
c
      character*16 sgmnam(nsa)
      real*8 x(nsa,3)
c
      call argos_diana_settag(sgmnam,byte_mb(i_tag),dbl_mb(i_val),x)
c
      return
      end
