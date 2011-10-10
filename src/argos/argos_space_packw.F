      subroutine argos_space_packw(numwm,iwl,iwlp)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer numwm,iwl(mwm,miw2),iwlp(mwm,npackw)
c
      integer i,j,k,ipck
c
      ipck=0
      do 1 j=1,npackw
      do 2 i=1,numwm
      iwlp(i,j)=0
    2 continue
      do 3 k=1,ipackw(j)
      ipck=ipck+1
      do 4 i=1,numwm
      iwlp(i,j)=lshift(iwlp(i,j),nbitw(ipck))+iwl(i,ipck)
    4 continue
    3 continue
    1 continue
c
      return
      end
c $Id$
