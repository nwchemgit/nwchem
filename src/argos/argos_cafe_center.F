      subroutine argos_cafe_center(xw,numw,isfrc,xs,nums,
     + idscb,nscb,icentr)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "msgids.fh"
c
      integer numw,nums,nscb,icentr
      integer isfrc(msa),idscb(nscb)
      real*8 xw(mwm,3,mwa),xs(msa,3)
c
      real*8 xb(4)
      integer i,isa,iwm,iwa,j,n,i0,i9
c
      if(nsa.le.0) return
c
      i0=1
      i9=3
      if(icentr.eq.1) i0=3
      if(icentr.eq.2) i9=2
c
      do 1 i=1,3
      xb(i)=zero
    1 continue
c
      n=0
      do 2 isa=1,nums
      do 3 j=1,nscb
      if(idscb(j).eq.isfrc(isa)) then
      n=n+1
      do 4 i=i0,i9
      xb(i)=xb(i)+xs(isa,i)
    4 continue
      endif
    3 continue
    2 continue
c
      xb(4)=dble(n)
c
      call ga_dgop(mcf_49,xb,4,'+')
c
      do 5 i=i0,i9
      xb(i)=xb(i)/xb(4)
    5 continue
c
      do 6 i=i0,i9
      do 7 isa=1,nums
      xs(isa,i)=xs(isa,i)-xb(i)
    7 continue
    6 continue
c
      do 8 i=i0,i9
      do 9 iwa=1,nwa
      do 10 iwm=1,numw
      xw(iwm,i,iwa)=xw(iwm,i,iwa)-xb(i)
   10 continue
    9 continue
    8 continue
c
      return
      end
c $Id$
