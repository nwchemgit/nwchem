      logical function pre_removw(xr,
     + xw,vw,mwm,mwa,nwm,nwa,xwc,vwc,mwmc,nwmc)
c
c $Id: pre_removw.F,v 1.1 2000-05-02 21:17:30 d3j191 Exp $
c
      implicit none
c
      real*8 xr(3,2)
      integer mwa,mwm,nwa,nwm,mwmc,nwmc
      real*8 xw(3,mwa,mwm),xwc(3,mwa,mwmc)
      real*8 vw(3,mwa,mwm),vwc(3,mwa,mwmc)
c
      integer i,j,k,l,new
      write(*,1000)
 1000 format(/,' Removing selected solvent')
c
      do 1 i=1,3
      if(xr(i,1).lt.xr(i,2)) then
      new=0
      do 2 j=1,nwm
      do 3 k=1,nwa
      if(xw(i,j,k).lt.xr(i,1)) goto 2
      if(xw(i,j,k).gt.xr(i,2)) goto 2
    3 continue
      new=new+1
      do 5 l=1,3
      do 6 k=1,nwa
      xw(l,k,new)=xw(l,k,i)
      vw(l,k,new)=vw(l,k,i)
    6 continue
    5 continue
    2 continue
      nwm=new
      new=0
      do 7 j=1,nwmc
      do 8 k=1,nwa
      if(xwc(i,j,k).lt.xr(i,1)) goto 7
      if(xwc(i,j,k).gt.xr(i,2)) goto 7
    8 continue
      new=new+1
      do 9 l=1,3
      do 10 k=1,nwa
      xwc(l,k,new)=xwc(l,k,i)
      vwc(l,k,new)=vwc(l,k,i)
   10 continue
    9 continue
    7 continue
      nwmc=new
      endif
    1 continue
c
      pre_removw=.true.
c
      return
      end
      
