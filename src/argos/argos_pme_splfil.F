      subroutine argos_pme_splfil(w,arr,darr)
c
      implicit none
c
#include "argos_pme_common.fh"
c
      real*8 arr(morder),darr(morder),w
c
      integer i,j
      real*8 rinv
c
      arr(morder)=zero
      arr(2)=w
      arr(1)=one-w
      do 1 i=3,morder-1
      rinv=one/dble(i-1)
      arr(i)=rinv*w*arr(i-1)
      do 2 j=1,i-2
      arr(i-j)=rinv*((w+dble(j))*arr(i-j-1)+(dble(i-j)-w)*arr(i-j))
    2 continue
      arr(1)=rinv*(one-w)*arr(1)
    1 continue
      darr(1)=-arr(1)
      do 3 j=2,morder
      darr(j)=arr(j-1)-arr(j)
    3 continue
      rinv=one/dble(morder-1)
      arr(morder)=rinv*w*arr(morder-1)
      do 4 j=1,morder-2
      arr(morder-j)=rinv*((w+dble(j))*arr(morder-j-1)+
     + (dble(morder-j)-w)*arr(morder-j))
    4 continue
      arr(1)=rinv*(one-w)*arr(1)
c
      return
      end
c $Id$
