      logical function pre_center(xw,mwm,mwa,nwm,nwa,xwc,mwmc,nwmc,
     + xs,msa,nsa)
c
c $Id: pre_center.F,v 1.7 1998-04-29 18:35:12 d3j191 Exp $
c
      implicit none
c
      integer mwa,mwm,nwa,nwm,msa,nsa,mwmc,nwmc
      real*8 xw(3,mwa,mwm),xwc(2,mwa,mwmc),xs(3,msa)
c
      integer i,j,k
      real*8 s
c
      do 1 k=1,3
      s=0.0d0
      do 2 i=1,nsa
      s=s+xs(k,i)
    2 continue
      s=s/dble(nsa)
      do 3 i=1,nsa
      xs(k,i)=xs(k,i)-s
    3 continue
      if(nwmc.gt.0) then
      do 4 i=1,nwmc
      do 5 j=1,nwa
      xwc(k,j,i)=xwc(k,j,i)-s
    5 continue
    4 continue
      endif
      if(nwm.gt.0) then
      do 6 i=1,nwm
      do 7 j=1,nwa
      xw(k,j,i)=xw(k,j,i)-s
    7 continue
    6 continue
      endif
    1 continue
c
      pre_center=.true.
      return
      end
      
