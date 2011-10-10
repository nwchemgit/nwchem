      subroutine argos_membrane_init(ismol,mm,xs,xsm,fm)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ismol(msa),mm(msa,2)
      real*8 xs(msa,3),xsm(msm,3),fm(msm,7)
c
      integer i,j,k,numg
      real*8 dx,dk
c
      numg=0
c
      do 1 i=1,msm
      mm(i,1)=0
    1 continue
      do 2 i=1,nsaloc
      mm(ismol(i),1)=mm(ismol(i),1)+1
      mm(i,2)=ismol(i)
    2 continue
      if(np.gt.1) call ga_dgop(mrg_d49,mm,msm,'+')
      do 3 i=1,nsaloc
      if(mm(ismol(i),1).eq.1) then
      k=0
      dk=zero
      do 4 j=1,nsm
      if(ismol(i).ne.j) then
      dx=(xs(i,1)-xsm(j,1))**2+(xs(i,2)-xsm(j,2))**2+
     + (xs(i,3)-xsm(j,3))**2
      if(k.eq.0.or.dx.lt.dk) then
      dk=dx
      k=j
      endif
      endif
    4 continue
      mm(i,2)=k
      numg=numg+1
      endif
    3 continue
      if(numg.gt.0.and.me.eq.0) then
      write(*,2000) numg
 2000 format(' Regrouping of',i5,' atoms',/)
      endif
c
      return
      end
c $Id$
