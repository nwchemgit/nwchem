      logical function pre_center(xw,mwm,mwa,nwm,nwa,xwc,mwmc,nwmc,
     + xs,msa,nsa)
c
c $Id: pre_center.F,v 1.10 2000-08-31 00:07:00 d3j191 Exp $
c
      implicit none
c
      integer mwa,mwm,nwa,nwm,msa,nsa,mwmc,nwmc
      real*8 xw(3,mwa,mwm),xwc(3,mwa,mwmc),xs(3,msa)
c
      integer i,j,k,imax,imin
      real*8 s,smin,smax
c
      write(*,1000)
 1000 format(/,' Centering solute')
c
      do 1 k=1,3
      smin=xs(k,1)
      smax=xs(k,1)
      imin=1
      imax=1
      do 2 i=2,nsa
      if(xs(k,i).lt.smin) imin=i
      if(xs(k,i).gt.smax) imax=i
      smin=min(smin,xs(k,i))
      smax=max(smax,xs(k,i))
    2 continue
      s=(smax+smin)/2.0d0
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
      write(*,'(i2,2i5,3f12.6)') k,imin,imax,smin,smax,s
    1 continue
c
      pre_center=.true.
      return
      end
      
