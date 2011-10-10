      subroutine argos_cafe_para(ia,nam,wgt,num)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer ia
      character*6 nam(24)
      real*8 wgt(24)
      integer num(24)
c
      if(ia.lt.1.or.ia.gt.mat) call md_abort('Error in para',0)
c
      call argos_cafe_copya(ia,nam,wgt,num,byte_mb(i_nam),
     + dbl_mb(i_mas),int_mb(i_num),int_mb(i_typ))
c
      if(lfree) then
      if(abs(wgt(1)-wgt(2)).gt.tiny) ip2(13)=.true.
      if(abs(wgt(1)-wgt(3)).gt.tiny) ip3(13)=.true.
      if(abs(wgt(2)-wgt(3)).gt.tiny) ith(13)=.true.
      endif
c
      return
      end
c
c $Id$
