      integer function argos_space_btop(ibox,ibownr)
c
      implicit none
c
#include "argos_space_common.fh"
c
      integer ibox,ibownr(maxbox,3)
      integer iboxx,iboxy,iboxz
c
      iboxx=mod(ibox,nbx)
      iboxy=mod((ibox-iboxx)/nbx,nby)
      iboxz=((ibox-iboxx)/nbx-iboxy)/nby
      argos_space_btop=(ibownr(iboxz+1,3)*npy+ibownr(iboxy+1,2))*npx+
     + ibownr(iboxx+1,1)
c
      return
      end
c $Id$
