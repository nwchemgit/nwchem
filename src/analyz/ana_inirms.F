      subroutine ana_inirms(xref,xrms,wt)
c
c $Id: ana_inirms.F,v 1.2 2000-03-16 20:45:37 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
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
