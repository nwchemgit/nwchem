      subroutine argos_cafe_parv(ia,ja,vdwi)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer ia,ja
      real*8 vdwi(24,4)
c
      integer i
c
c     van der Waals parameters are entered here as:
c
c     vdwai(1:6) : c6a (set 1,2,3),  c6b (set 1,2,3)
c     vdwri(1:6) : c12a (set 1,2,3), c12b (set 1,2,3)
c
      if(ia.lt.1.or.ia.gt.mat) call md_abort('Error in para',0)
      if(ja.lt.1.or.ja.gt.mat) call md_abort('Error in para',0)
c
      call argos_cafe_copyv(ia,ja,vdwi,dbl_mb(i_vdw))
c
      if(lfree) then
      do 1 i=1,4
      if(abs(vdwi(1,i)-vdwi(2,i)).gt.tiny) ip2(14)=.true.
      if(abs(vdwi(1,i)-vdwi(3,i)).gt.tiny) ip3(14)=.true.
      if(abs(vdwi(2,i)-vdwi(3,i)).gt.tiny) ith(14)=.true.
    1 continue
      endif
c
      return
      end
c $Id$
