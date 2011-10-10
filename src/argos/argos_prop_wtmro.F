      subroutine argos_prop_wtmro(lfnmro,ndec,dfr,dfrm)
c
      implicit none
c
#include "argos_prop.fh"
c
      integer lfnmro,ndec
      real*8 dfr(mda),dfrm(mda)
c
      integer i
c
      if(me.eq.0) then
      write(lfnmro) nderiv,nda,nprop,nsum,nsumt,ndec,maxp
      write(lfnmro) deriv
      write(lfnmro) (dfr(i),i=1,nda)
      write(lfnmro) (dfrm(i),i=1,nda)
      write(lfnmro) tsum,t2sum,tsumt,t2sumt
      write(lfnmro) (psum(i),i=1,maxp)
      write(lfnmro) (p2sum(i),i=1,maxp)
      write(lfnmro) (pslop(i),i=1,maxp)
      write(lfnmro) (psumt(i),i=1,maxp)
      write(lfnmro) (p2sumt(i),i=1,maxp)
      write(lfnmro) (pslopt(i),i=1,maxp)
      endif
      if(ndec.gt.0) call argos_cafe_wrtmro(lfnmro)
c
      return
      end
c $Id$
