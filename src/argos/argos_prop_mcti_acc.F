      logical function argos_prop_mcti_acc(ida,dfr,dfrm,fdata,
     + aver,drift,stderr,corerr,ratio)
c
      implicit none
c
#include "argos_prop.fh"
c
      integer ida
      real*8 dfr(mda),dfrm(mda),fdata(28)
      real*8 aver,drift,stderr,corerr,ratio,cerror
c
      integer i
      logical done
      real*8 dfrnom
c
      dfree=zero
      dfrnom=zero
      do 1 i=1,24
      dfree=dfree+fdata(i)
      if(i.ne.1.and.i.ne.13) dfrnom=dfrnom+fdata(i)
      deriv(i)=deriv(i)+fdata(i)
    1 continue
      nderiv=nderiv+1
c
      dfr(ida)=dfree
      dfrm(ida)=dfrnom
      if(dfree.ne.dfrnom) lfreem=.true.
      nda=ida
c
    2 continue
c
      if(lerror) then
      call error(lauto,lappr,1000,dfr,ida,
     + aver,drift,stderr,corerr,ratio)
      cerror=corerr
      if(.not.lauto) cerror=samrat*corerr
      if(.not.lauto.and.cerror.lt.edacq) then
      lauto=.true.
      lappr=.true.
      goto 2
      else
      if(lauto) samrat=ratio
      lauto=.false.
      lappr=.false.
      endif
      done=cerror.lt.edacq.and.drift.lt.ddacq
      else
      done=.false.
      endif
c
      argos_prop_mcti_acc=done
      return
      end
c $Id$
