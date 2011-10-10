      subroutine argos_cafe_slowt(vw,numw,vs,nums,iwat,isat,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "msgids.fh"
c
      integer numw,nums
      real*8 vw(mwm,3,mwa),vs(msa,3),wgt(mat,mset)
      integer iwat(mwa),isat(msa)
c
      integer i,j,k
      real*8 vb(3)
c
      do 1 i=1,3
      vb(i)=zero
    1 continue
c
      do 2 k=1,3
      do 3 j=1,nwa
      do 4 i=1,numw
      vb(k)=vb(k)+wgt(iwat(j),iset)*vw(i,k,j)
    4 continue
    3 continue
      do 5 i=1,nums
      vb(k)=vb(k)+wgt(isat(i),iset)*vs(i,k)
    5 continue
      vb(k)=vb(k)/wbox
    2 continue
c
      call ga_dgop(mcf_50,vb,3,'+')
c
      do 6 k=1,3
      do 7 j=1,nwa
      do 8 i=1,numw
      vw(i,k,j)=vw(i,k,j)-vb(k)
    8 continue
    7 continue
      do 9 i=1,nums
      vs(i,k)=vs(i,k)-vb(k)
    9 continue
    6 continue
c
      return
      end
c $Id$
