      subroutine argos_space_unpackw(numwm,iwl,iwlp)
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
      ipck=miw2+1
      do 1 j=npackw,1,-1
      do 2 k=ipackw(j),1,-1
      ipck=ipck-1
      do 3 i=1,numwm
      iwl(i,ipck)=iand(iwlp(i,j),miwl(ipck))
      iwlp(i,j)=ishft(iwlp(i,j)-iwl(i,ipck),-nbitw(ipck))
    3 continue
    2 continue
    1 continue
c
      return
      end
c $Id$
