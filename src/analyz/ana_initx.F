      subroutine ana_initx(nx,ix,x)
c
c $Id: ana_initx.F,v 1.1 1998-12-29 18:01:03 d3j191 Exp $
c
c     reference coordinated entry of analysis module:
c
c     integer nx         : number of atoms
c     integer mx         : leading dimension atomic arrays
c     integer ix(1:nx)   : global atom numbers
c     real*8 x(1:nx,1:3) : atomic coordinates
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      real*8 t
      integer nx
      integer ix(msa)
      real*8 x(msa,3)
c
      call ana_cpxref(ix,x,nx,dbl_mb(i_xref))
c
      return
      end
      subroutine ana_cpxref(ix,x,nx,xref)
c
      implicit none
c
#include "ana_params.fh"
#include "ana_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nx
      integer ix(nx)
      real*8 x(msa,3),xref(nsa,3)
c
      integer i,j
c
      do 1 j=1,3
      do 2 i=1,nsa
      xref(i,j)=zero
    2 continue
    1 continue
c
      do 3 i=1,nx
      do 4 j=1,3
      xref(ix(i),j)=x(i,j)
    4 continue
    3 continue
c
      if(np.gt.1) call ga_dgop(mag_d01,xref,3*nsa,'+')
c
      return
      end
