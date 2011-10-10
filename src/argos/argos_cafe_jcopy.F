      subroutine argos_cafe_jcopy(idx,jdx,n,ip,nip)
c
      implicit none
c
      integer n,nip
      integer idx(n),jdx(n),ip(0:nip,2)
c
      integer i
c
      do 1 i=1,n
      ip(i,1)=-idx(i)
      ip(i,2)=jdx(i)
    1 continue
      if(n.lt.nip) then
      do 2 i=n+1,nip
      ip(i,1)=0
    2 continue
      endif
c
      return
      end
c $Id$
