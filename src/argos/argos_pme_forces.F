      subroutine argos_pme_forces(fw,nwmloc,fs,isq,isgm,nsaloc)
c
      implicit none
c
#include "argos_pme_common.fh"
#include "mafdecls.fh"
c
      integer nwmloc,nsaloc
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer isq(msa),isgm(msa)
c
      call argos_pme_for(fw,nwmloc,fs,isq,isgm,nsaloc,int_mb(i_lnod),
     + int_mb(i_kfr),int_mb(i_kto))
c
      return
      end
c $Id$
