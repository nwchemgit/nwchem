      subroutine argos_zinit(iwz,isz)
c
      implicit none
c
#include "argos.fh"
c
      integer iwz(mwm),isz(msa)
c
      integer i
c
      do 1 i=1,mwm
      iwz(i)=0
    1 continue
c
      do 2 i=1,msa
      isz(i)=0
    2 continue
c
      return
      end
c $Id$
