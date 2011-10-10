      subroutine argos_cafe_copwiq(id,ia,iq,iwa,iwq,wgt,chg,vdw)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer i,j,id,ia,iq,iwa(mwa),iwq(mwa)
      real*8 wgt(mat,mset),chg(mqt,mqp,mset),vdw(mat,mat,map,mset)
c
      iwa(id)=ia
      iwq(id)=iq
c
      if(lfree) then
      if(abs(wgt(ia,1)-wgt(ia,2)).gt.tiny) then
      ip2(1)=.true.
      ith(1)=.true.
      endif
      do 1 j=1,map
      do 2 i=1,mat
      if(abs(vdw(ia,i,j,1)-vdw(ia,i,j,2)).gt.tiny) then
      ip2(2)=.true.
      ith(2)=.true.
      endif
      if(nparms.gt.2) then
      if(abs(vdw(ia,i,j,1)-vdw(ia,i,j,3)).gt.tiny) ip3(2)=.true.
      endif
    2 continue
    1 continue
      if(abs(chg(iq,1,1)-chg(iq,1,2)).gt.tiny) then
      ip2(4)=.true.
      ith(4)=.true.
      endif
      if(abs(chg(iq,2,1)-chg(iq,2,2)).gt.tiny) then
      ip2(5)=.true.
      ith(5)=.true.
      endif
      if(nparms.gt.2) then
      if(abs(wgt(ia,1)-wgt(ia,3)).gt.tiny) ip3(1)=.true.
      if(abs(chg(iq,1,1)-chg(iq,1,3)).gt.tiny) ip3(4)=.true.
      if(abs(chg(iq,2,1)-chg(iq,2,3)).gt.tiny) ip3(5)=.true.
      endif
      endif
c
      return
      end
c $Id$
