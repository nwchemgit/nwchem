      logical function argos_prop_mcti_step(ida,lda)
c
      implicit none
c
#include "argos_prop.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      logical argos_prop_mcti_acc
      external argos_prop_mcti_acc
c
      integer ida,lda
c
      real*8 fdata(28)
      real*8 aver,drift,stderr,corerr,ratio
      logical done
c
      lerror=ida.gt.lda
c
      call argos_cafe_mcti(fdata)
c
      call ga_dgop(mrg_d44,fdata,28,'+')
c
      if(me.eq.0) then
      done=argos_prop_mcti_acc(ida,dbl_mb(i_dfr),dbl_mb(i_dfrm),fdata,
     + aver,drift,stderr,corerr,ratio)
      endif
c
      if(np.gt.1) then
      call ga_brdcst(mrg_d46,done,ma_sizeof(mt_log,1,mt_byte),0)
      endif
c
      argos_prop_mcti_step=done
c
      return
      end
c $Id$
