      subroutine argos_pme_wait(iflag,icolor)
c
      implicit none
c
#include "argos_pme_common.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 timer_wall
      external timer_wall
c
      integer iflag,icolor
      integer il,ih,jl,jh,icolr
c
      call ga_distribution(ga_flg,me,il,ih,jl,jh)
    1 continue
      call ga_get(ga_flg,il,ih,jl+iflag,jl+iflag,icolr,1)
      if(icolr.ne.icolor) goto 1
c
      return
      end
c $Id$
