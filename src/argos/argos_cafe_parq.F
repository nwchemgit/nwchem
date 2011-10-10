      subroutine argos_cafe_parq(ia,chgi)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer ia
      real*8 chgi(24,3)
c
      integer i
c
      if(ia.lt.1.or.ia.gt.mqt) call md_abort('Error in parq',0)
c
      call argos_cafe_copyq(ia,chgi,dbl_mb(i_chg))
c
      if(lfree) then
      do 1 i=1,2
      if(abs(chgi(1,i)-chgi(2,i)).gt.tiny) ip2(15+i)=.true.
      if(abs(chgi(1,i)-chgi(3,i)).gt.tiny) ip3(15+i)=.true.
      if(abs(chgi(2,i)-chgi(3,i)).gt.tiny) ith(15+i)=.true.
    1 continue
      endif
c
      return
      end
c $Id$
