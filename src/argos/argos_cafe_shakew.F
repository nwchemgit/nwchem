      subroutine argos_cafe_shakew(xw,yw,iwgm,iwdt,numw)
c $Id$
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer numw
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      integer iwgm(mwm),iwdt(mwm)
c
      call argos_cafe_shkw(numw,iwgm,iwdt,xw,yw,
     + mbt(1),numb(1),mbp(1),int_mb(i_ibnd(1)),dbl_mb(i_bnd(1)),
     + int_mb(i_iwa),dbl_mb(i_mas),
     + dbl_mb(i_s3r1),dbl_mb(i_s1r1),dbl_mb(i_s1r2),dbl_mb(i_s1r3),
     + dbl_mb(i_s1r4),dbl_mb(i_smr1))
c
      return
      end
