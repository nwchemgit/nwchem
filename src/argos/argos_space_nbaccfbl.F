      subroutine argos_space_nbaccfbl(ibbl,lhandl,fw,fs,lpair,iwz,isz)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
c
      integer ibbl,lhandl
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      logical lpair
      integer iwz(mwm),isz(msa)
c
      call argos_space_nbacfbl(ibbl,lhandl,fw,fs,int_mb(i_bb),
     + int_mb(i_jpl),lpair,iwz,isz)
c
      return
      end
c $Id$
