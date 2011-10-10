      subroutine argos_membrane_forces(mm,fm,xs,xsm,fs,ws)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer mm(msa,2)
      real*8 fm(msm,7),xs(msa,3),xsm(msm,3),fs(msa,3),ws(msa)
c
      integer i,j
      real*8 factor
c
      do 1 i=1,msm
      fm(i,1)=zero
      fm(i,2)=zero
      fm(i,3)=zero
      fm(i,4)=zero
      fm(i,5)=zero
      fm(i,6)=zero
      fm(i,7)=zero
    1 continue
c
      do 2 i=1,nsaloc
      factor=one/ws(i)
      fm(mm(i,2),1)=fm(mm(i,2),1)+factor*fs(i,1)
      fm(mm(i,2),2)=fm(mm(i,2),2)+factor*fs(i,2)
      fm(mm(i,2),3)=fm(mm(i,2),3)+factor*fs(i,3)
      fm(mm(i,2),4)=fm(mm(i,2),4)+factor*
     + ((xs(i,1)-xsm(mm(i,2),1))*fs(i,2)-
     +  (xs(i,2)-xsm(mm(i,2),2))*fs(i,1))
    2 continue
      if(np.gt.1) call ga_dgop(mrg_d50,fm,4*msm,'+')
c
      do 3 i=1,nsm
      fm(i,4)=fm(i,3)
      if(me.eq.0) write(lfnout,1000) i,(fm(i,j)/dble(mm(i,1)),j=1,3)
 1000 format(i5,3f12.3)
    3 continue
c
c     molecular rotations only
c
      if(imembr.eq.2) then
      do 4 i=1,nsaloc
      fs(i,1)=fs(i,1)-ws(i)*fm(mm(i,2),1)/dble(mm(i,1))
      fs(i,2)=fs(i,2)-ws(i)*fm(mm(i,2),2)/dble(mm(i,1))
      fs(i,3)=fs(i,3)-ws(i)*fm(mm(i,2),3)/dble(mm(i,1))
    4 continue
      endif
c
      return
      end
c $Id$
