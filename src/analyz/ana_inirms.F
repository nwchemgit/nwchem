      subroutine ana_inirms(xref,xrms,wt)
c
c $Id: ana_inirms.F,v 1.1 1998-12-29 18:01:02 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
c
      real*8 xref(nsa,3),xrms(nsa),wt(nsa)
c
      integer i,j
c
      do 1 i=1,nsa
      xrms(i)=zero
      wt(i)=one
    1 continue
c
      return
      end
