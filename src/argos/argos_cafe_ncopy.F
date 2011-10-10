      subroutine argos_cafe_ncopy(idx,jdx,n,ip,nip)
c
      implicit none
c
      integer n,nip
      integer idx(n),jdx(n),ip(0:nip,2)
c
      integer i
c
      do 1 i=0,nip
      ip(i,1)=0
      ip(i,2)=0
    1 continue
      do 2 i=1,n
      ip(idx(i),1)=i
      ip(i,2)=jdx(i)
    2 continue
      do 3 i=1,nip
      if(ip(i,1).eq.0) ip(i,1)=ip(i-1,1)
    3 continue
c
      return
      end
c $Id$
