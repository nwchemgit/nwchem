      subroutine argos_diana_inirms(xref,xrms,wt)
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
c
      real*8 xref(nsa,3),xrms(nsa),wt(nsa)
c
      integer i
c
      do 1 i=1,nsa
      xrms(i)=zero
      wt(i)=one
    1 continue
c
      return
      end
