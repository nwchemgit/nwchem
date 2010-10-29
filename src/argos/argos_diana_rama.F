      subroutine argos_diana_rama()
c
c $Id$
c
      implicit none
c
#include "argos_diana_params.fh"
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      character*255 string
      integer ilt,iht,jlt,jht,ndxt,ldt
c
      string=card(8:80)
      strram=string
c
c     get access to local trajectory and time
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
      call ga_distribution(ga_trt,me,ilt,iht,jlt,jht)
      call ga_access(ga_trt,ilt,iht,jlt,jht,ndxt,ldt)
c
      call argos_diana_rama1(int_mb(i_iram),dbl_mb(ndx),dbl_mb(ndxt))
c
      return
      end
