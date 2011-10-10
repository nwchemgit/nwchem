      subroutine argos_cafe_fnrm(wgt,iwat,iwdt,fw,numw,
     + isat,isdt,fs,nums,fnorm,fmax)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      real*8 wgt(mat,mset)
      real*8 fw(mwm,3,mwa),fs(msa,3),fnorm,fmax
      integer iwat(mwa),isat(msa),numw,nums
      integer isdt(msa),iwdt(mwm)
c
      integer i,j
      real*8 fm,factor
c
      fnorm=zero
      fmax=zero
c
      if(numw.gt.0) then
      do 1 j=1,nwa
      factor=(one/wgt(iwat(j),iset))*(one/wgt(iwat(j),iset))
      do 2 i=1,numw
      if(iand(iwdt(i),lfixed).eq.0) then
      fm=fw(i,1,j)**2+fw(i,2,j)**2+fw(i,3,j)**2
      fnorm=fnorm+fm
      fmax=max(fmax,factor*fm)
      endif
    2 continue
    1 continue
      endif
c
      if(nums.gt.0) then
      do 3 i=1,nums
      if(iand(isdt(i),lfixed).eq.0) then
      fm=fs(i,1)**2+fs(i,2)**2+fs(i,3)**2
      fnorm=fnorm+fm
      factor=(one/wgt(isat(i),iset))*(one/wgt(isat(i),iset))
      fmax=max(fmax,factor*fm)
      endif
    3 continue
      endif
c
      if(np.gt.1) then
      call ga_dgop(mcf_44,fnorm,1,'+')
      call ga_dgop(mcf_45,fmax,1,'max')
      endif
c
      fnorm=sqrt(fnorm)
      fmax=sqrt(fmax)
c
      return
      end
c $Id$
