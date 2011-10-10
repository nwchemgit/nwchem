      subroutine argos_cafe_copyq(ia,chgi,chg)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer ia
      real*8 chgi(24,3),chg(mqt,mqp,mset)
c
      integer i
c
      if(.not.lfree) then
      do 1 i=1,nparms
      chg(ia,1,i)=qfac*chgi(i,1)
      chg(ia,2,i)=chgi(i,2)
      chg(ia,3,i)=qfac*chgi(i,3)
    1 continue
      else
      do 2 i=1,3
      chg(ia,1,i)=qfac*chgi(i,1)
      chg(ia,2,i)=chgi(i,2)
      chg(ia,3,i)=qfac*chgi(i,3)
    2 continue
      do 3 i=1,2
      chg(ia,i,4)=chg(ia,i,3)-chg(ia,i,2)
      chg(ia,i,5)=chg(ia,i,2)
      chg(ia,i,6)=chg(ia,i,3)
    3 continue
      endif
c
      return
      end
c $Id$
