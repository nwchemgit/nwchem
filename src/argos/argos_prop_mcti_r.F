      subroutine argos_prop_mcti_r(rlambd,dlambd,ndec,dfr,dfrm,taver,
     + ep2ave,ep3ave,ep2avm,ep3avm,ebias,dfbias)
c
      implicit none
c
#include "argos_prop.fh"
c
      real*8 rlambd,dlambd,dfr(mda),dfrm(mda)
      real*8 taver,ep2ave,ep3ave,ebias,dfbias,ep2avm,ep3avm
      integer i,ndec
c
      if(me.eq.0) then
      write(lfngib,1000) nderiv,nda,rlambd,dlambd,ndec,nsa,
     + nbias,ebias/dble(nsumt),lfreem
 1000 format(2i7,2f12.6,2i8,i4,e20.12,4x,l1)
      write(lfngib,1001) deriv
 1001 format(4e20.12)
      write(lfngib,1002) (dfr(i),i=1,nda)
      write(lfngib,1002) (dfrm(i),i=1,nda)
 1002 format(4e20.12)
      write(lfngib,1003) nsumt,taver/dble(nsumt),
     + ep2ave/dble(nsumt),ep3ave/dble(nsumt),dfbias/ebias,
     + ep2avm/dble(nsumt),ep3avm/dble(nsumt)
 1003 format(i10,/,4e20.12,/,2e20.12)
      endif
      if(ndec.gt.0) call argos_cafe_wrtgib(lfngib)
c
      return
      end
c $Id$
