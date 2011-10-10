      subroutine argos_diana_supimp(xtrj)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "argos_diana_common.fh"
c
      real*8 zero
      parameter(zero=0.0d0)
c
      integer nfr,nxa
      real*8 xtrj(ntrj,3,nfrdim)
c
      integer i,j,k
c
      nwm=0
      do 1 k=1,nfrme
      call super(xtrj(1,1,k),ntrj,ntrj,
     + dbl_mb(i_xref),dbl_mb(i_swt),dbl_mb(i_wt),dbl_mb(i_xrms),
     + nsa,msa,.true.,rms0,rms1,dbl_mb(i_wdat),mwm,nwm,mwa,nwa,.false.)
    1 continue
c
      return
      end
c $Id$
