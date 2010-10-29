      subroutine ana_rmsdev(x,ix,nx,xref,wt,wts,xrms)
c $Id$
      implicit none
c
#include "ana_params.fh"
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer nx
      integer ix(nx)
      real*8 x(msa,3),xref(nsa,3),wt(nsa),wts(nsa),xrms(nsa)
c
      real*8 rms0,rms1
c
      call super2(x,ix,nx,msa,xref,wt,wts,xrms,nsa,nsa,.false.,
     + rms0,rms1,dbl_mb(i_wdat),mwm,nwm,mwa,nwa,.true.)
c
      if(me.eq.0) then
      write(lfnrms,1000) time,rms0,rms1
 1000 format(3f12.6)
      endif
c
      return
      end
