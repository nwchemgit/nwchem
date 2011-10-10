      subroutine argos_diana_histo(x,w,ihi)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ihi,ito
      real*8 x(nsa,3),w(mwm,mwa,3)
c
      integer i,ia,j,ihndx
c
      dhis=1.2d0*box(3)/real(idhis(ihi,2))
c
      ito=ldef(idhis(ihi,1))
      if(ito.gt.0) then
      do 1 i=1,ito
      ia=idef(idhis(ihi,1),i)
      ihndx=(x(ia,3)-rhis)/dhis
      ihis(ihndx,ihi)=ihis(ihndx,ihi)+1
    1 continue
      else
      do 2 j=1,nwm
      do 3 i=1,-ito
      ia=idef(idhis(ihi,1),i)
      ihndx=(w(j,ia,3)-rhis)/dhis
      ihis(ihndx,ihi)=ihis(ihndx,ihi)+1
    3 continue
    2 continue
      endif
c
      return
      end
c $Id$
