      subroutine argos_cafe_mcti_w(isat,isga,vs,nums,wgt,dera)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer nums
      integer isat(msa),isga(msa)
      real*8 vs(msa,3),wgt(mat,mset)
      real*8 dera(6,nsatot)
c
      integer i
      real*8 vsq
c
      do 1 i=1,nums
      vsq=half*(vs(i,1)*vs(i,1)+vs(i,2)*vs(i,2)+vs(i,3)*vs(i,3))
      deriv(13,1)=deriv(13,1)-wgt(isat(i),4)*vsq
      if(npgdec.gt.1) dera(6,isga(i))=dera(6,isga(i))-wgt(isat(i),4)*vsq
      if(ipert2.and.ip2(13)) ep2(3)=ep2(3)+
     + vsq*wgt(isat(i),1)*(wgt(isat(i),1)/wgt(isat(i),2)-one)
      if(ipert3.and.ip3(13)) ep3(3)=ep3(3)+
     + vsq*wgt(isat(i),1)*(wgt(isat(i),1)/wgt(isat(i),3)-one)
    1 continue
c
      return
      end
c $Id$
