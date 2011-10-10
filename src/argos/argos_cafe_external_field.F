      subroutine argos_cafe_external_field(nums,fs,isq,numw,fw,iwq,chg)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer nums,numw
      integer isq(msa),iwq(mwa)
      real*8 fs(msa,3),fw(mwm,3,mwa),chg(mqt,mqp,mset)
c
      integer i,j,k
      real*8 efield,q,df
c
      efield=field*cos(2.0d-6*pi*ffreq*stime)
c
      if(numw.gt.0) then
      do 1 k=1,nwa
      q=chg(iwq(k),1,iset)*efield
      do 2 j=1,3
      df=q*fvect(j)
      do 3 i=1,numw
      fw(i,j,k)=fw(i,j,k)+df
    3 continue
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      do 4 j=1,3
      do 5 i=1,nums
      fs(i,j)=fs(i,j)+chg(isq(i),1,iset)*efield*fvect(j)
    5 continue
    4 continue
      endif
c
      return
      end
c $Id$
