      subroutine argos_diana_rdf(isel,xs,xw,irdf)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer isel(msa),irdf(nsel,mwa,nrdf)
      real*8 xs(msa,3),xw(mwm,mwa,3)
      integer i,j,k,l,m
      real*8 d
c
      i=0
      do 1 j=1,nsa
      if(isel(j).gt.0) then
      i=i+1
      do 2 k=1,nwm
      do 3 l=1,nwa
      d=sqrt((xs(j,1)-xw(k,l,1))**2+(xs(j,2)-xw(k,l,2))**2+
     + (xs(j,3)-xw(k,l,3))**2)
      m=int(dble(nrdf*rrdf)/d)
      if(m.le.nrdf) irdf(i,l,m)=irdf(i,l,m)+1
    3 continue
    2 continue
      endif
    1 continue
c
      numrdf=numrdf+1
c
      return
      end
c $Id$
