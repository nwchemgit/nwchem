      subroutine dia_initx
c
c $Id$
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
#include "dia_common.fh"
#include "mafdecls.fh"
c
c
c
      call dia_rdref(dbl_mb(i_xref))
c
      return
      end

