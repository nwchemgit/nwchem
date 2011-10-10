      logical function argos_cafe_shakep(xs,ys,igan,isat,isdt,ishp,nums)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nums
      real*8 xs(msa,3),ys(msa,3)
      integer igan(msa),isat(msa),isdt(msa),ishp(msa)
      integer nitall
c
      argos_cafe_shakep=.true.
c
      if(lpmf) then
      call argos_cafe_shkp(nums,xs,ys,igan,isat,ishp,
     + int_mb(i_ipmf),int_mb(i_jpmf),
     + dbl_mb(i_rpmf),dbl_mb(i_xpmf),dbl_mb(i_ypmf),dbl_mb(i_mas))
      nitall=niters-1
      call ga_igop(mcf_64,nitall,1,'+')
      argos_cafe_shakep=nitall.le.0
      endif
c
      return
      end
c $Id$
