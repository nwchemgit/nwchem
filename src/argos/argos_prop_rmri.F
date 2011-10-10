      logical function argos_prop_rmri(lfnmri,ndec,mropt,dfr,dfrm)
c
      implicit none
c
#include "argos_prop.fh"
c
      logical argos_cafe_rdmri
      external argos_cafe_rdmri
c
      integer lfnmri,ndec,mropt
      real*8 dfr(mda),dfrm(mda),dtmp
c
      integer i,nprp,nmxp
c
      read(lfnmri,err=9,end=9) nderiv,nda,nprp,nsum,nsumt,ndec,nmxp
      if(mropt.ne.2) then
      if(nprop.ne.nprp) then
      call md_abort('Number of properties nprop changed',0)
      endif
      if(maxp.ne.nmxp) then
      call md_abort('Number of properties maxp changed',0)
      endif
      endif
      read(lfnmri,err=9,end=9) deriv
      if(mda.lt.nda) then
      read(lfnmri,err=9,end=9) (dtmp,i=1,nda)
      read(lfnmri,err=9,end=9) (dtmp,i=1,nda)
      else
      read(lfnmri,err=9,end=9) (dfr(i),i=1,nda)
      read(lfnmri,err=9,end=9) (dfrm(i),i=1,nda)
      endif
      read(lfnmri,err=9,end=9) tsum,t2sum,tsumt,t2sumt
      if(mropt.ne.2) then
      read(lfnmri,err=9,end=9) (psum(i),i=1,maxp)
      read(lfnmri,err=9,end=9) (p2sum(i),i=1,maxp)
      read(lfnmri,err=9,end=9) (pslop(i),i=1,maxp)
      read(lfnmri,err=9,end=9) (psumt(i),i=1,maxp)
      read(lfnmri,err=9,end=9) (p2sumt(i),i=1,maxp)
      read(lfnmri,err=9,end=9) (pslopt(i),i=1,maxp)
      else
      do 1 i=1,6
      read(lfnmri,err=9,end=9) dtmp
    1 continue
      endif
      if(ndec.gt.0) then
      argos_prop_rmri=argos_cafe_rdmri(lfnmri)
      else
      argos_prop_rmri=.true.
      endif
      return
c
    9 continue
      argos_prop_rmri=.false.
      return
      end
c $Id$
