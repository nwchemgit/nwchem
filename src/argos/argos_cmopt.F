      subroutine argos_cmopt(vs,fs,fcm,ismol,isatt)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 vs(msa,3),fs(msa,3,2),fcm(msm,5)
      integer ismol(msa),isatt(msa)
c
      call argos_cafe_scmfor(icmopt,ismol,isatt,vs,fs,nsaloc,fcm)
c
      return
      end
c $Id$
