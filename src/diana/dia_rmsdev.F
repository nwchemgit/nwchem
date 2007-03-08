      subroutine dia_rmsd
c
c $Id: dia_rmsdev.F,v 1.5 2007-03-08 23:07:21 d3j191 Exp $
      implicit none
c
#include "dia_params.fh"
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
c
      if(card(8:13).eq.'averag') then
      elseif(card(8:13).eq.'refere') then
      call dia_rmsdev(dbl_mb(ndx),nsa,
     + min((me+1)*nfrme,nfrtot)-me*nfrdim,
     + dbl_mb(i_xref),dbl_mb(i_swt),dbl_mb(i_wt),dbl_mb(i_xrms))
      else
      call md_abort('Error in rmsdev directive',0)
      endif
c
      return
      end
      subroutine dia_rmsdev(xtrj,ntrj,nfr,xr,wt,wts,xrms)
c $Id: dia_rmsdev.F,v 1.5 2007-03-08 23:07:21 d3j191 Exp $
      implicit none
c
#include "dia_params.fh"
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ntrj,nfr
      real*8 xtrj(ntrj,3,nfr)
      real*8 xr(nsa,3),wt(nsa),wts(nsa),xrms(nsa)
c
      integer k
c
      do 1 k=1,nfr
      call super(xtrj(1,1,k),ntrj,ntrj,xr,wt,wts,xrms,nsa,nsa,
     + .false.,rms0,rms1,dbl_mb(i_wdat),mwm,nwm,mwa,nwa,.false.)
c
      write(*,1000) me,rms0,rms1
 1000 format(i5,2f12.6)
c
    1 continue
c
      return
      end
