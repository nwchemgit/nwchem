      subroutine argos_space_accpbl(ibbl,lhandl,pw,pwp,ps,psp,
     + lpair,iwz,isz)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      call argos_space_acpbl(ibbl,lhandl,pw,pwp,ps,psp,int_mb(i_bb),
     + int_mb(i_jpl))
c
      return
      end
c $Id$
