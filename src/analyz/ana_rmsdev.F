      subroutine ana_rmsdev(x,ix,nx,xref,wt,xrms)
c $Id: ana_rmsdev.F,v 1.6 2000-08-03 00:25:50 d3j191 Exp $
      implicit none
c
#include "ana_params.fh"
#include "ana_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nx
      integer ix(nx)
      real*8 x(msa,3),xref(nsa,3),wt(nsa),xrms(nsa)
c
      real*8 rms0,rms1
c
      call super(x,ix,nx,msa,xref,wt,xrms,nsa,nsa,.false.,rms0,rms1,
     + dbl_mb(i_wdat),mwm,nwm,mwa,nwa)
c
      if(me.eq.0) then
      write(lfnrms,1000) time,rms0,rms1
 1000 format(3f12.6)
      endif
c
      return
      end
