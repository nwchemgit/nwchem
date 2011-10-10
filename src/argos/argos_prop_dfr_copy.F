      integer function argos_prop_dfr_copy(ncopy,dfr,dfrm)
c
      implicit none
c
#include "argos_prop.fh"
c
      integer ncopy
      real*8 dfr(mda),dfrm(mda)
c
      integer i,ioff
c
      if(ncopy.lt.nda) then
      ioff=nda-ncopy
      do 1 i=1,ncopy
      dfr(i)=dfr(i+ioff)
      dfrm(i)=dfrm(i+ioff)
    1 continue
      nda=ncopy
      endif
c
      argos_prop_dfr_copy=nda
      return
      end
c $Id$
