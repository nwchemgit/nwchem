      subroutine dia_data(t,nx,ix,x)
c
c $Id$
c
c     data entry of analysis module:
c
c     real*8 t           : simulation time in ps
c     integer nx         : number of atoms
c     integer mx         : leading dimension atomic arrays
c     integer ix(1:nx)   : global atom numbers
c     real*8 x(1:nx,1:3) : atomic coordinates
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
c
      real*8 t
      integer nx
      integer ix(msa)
      real*8 x(msa,3)
c
      time=t
      ndata=ndata+1
c
      if(irmsd.gt.0) then
      call dia_rmsdev(x,ix,nx,dbl_mb(i_xref),dbl_mb(i_wt),
     + dbl_mb(i_xrms))
      endif
c
      return
      end
