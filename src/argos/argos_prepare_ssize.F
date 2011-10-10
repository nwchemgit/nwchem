      logical function argos_prepare_ssize(isgm,xs,msa,nsa,rsgm)
c
      implicit none
c
      integer msa,nsa
      integer isgm(msa)
      real*8 xs(3,msa)
      real*8 rsgm
c
      integer i,j
      real*8 d
c
      rsgm=0.0d0
c
      do 1 i=1,nsa-1
      do 2 j=i+1,nsa
      if(isgm(i).eq.isgm(j)) then
      d=(xs(1,i)-xs(1,j))**2+(xs(2,i)-xs(2,j))**2+(xs(3,i)-xs(3,j))**2
      rsgm=max(rsgm,d)
      endif
    2 continue
    1 continue
      if(rsgm.gt.0) rsgm=dsqrt(rsgm)
c
      argos_prepare_ssize=.true.
c
      return
      end
c $Id$
