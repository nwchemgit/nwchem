      logical function pre_removw(ir,xr,
     + xw,vw,mwm,mwa,nwm,nwa,xwc,vwc,mwmc,nwmc)
c
c $Id: pre_removw.F,v 1.2 2000-05-02 23:01:01 d3j191 Exp $
c
      implicit none
c
      integer ir
      real*8 xr(3,2)
      integer mwa,mwm,nwa,nwm,mwmc,nwmc
      real*8 xw(3,mwa,mwm),xwc(3,mwa,mwmc)
      real*8 vw(3,mwa,mwm),vwc(3,mwa,mwmc)
c
      integer i,j,k,l,idel,new
      write(*,1000)
 1000 format(/,' Removing selected solvent')
c
      if(ir.eq.0) then
      do 1 i=1,3
      if(xr(i,1).lt.xr(i,2)) then
      new=0
      do 2 k=1,nwm
      do 3 j=1,nwa
      if(xw(i,j,k).lt.xr(i,1)) goto 2
      if(xw(i,j,k).gt.xr(i,2)) goto 2
    3 continue
      new=new+1
      do 4 l=1,3
      do 5 j=1,nwa
      xw(l,j,new)=xw(l,j,k)
      vw(l,j,new)=vw(l,j,k)
    5 continue
    4 continue
    2 continue
      nwm=new
      new=0
      do 6 k=1,nwmc
      do 7 j=1,nwa
      if(xwc(i,j,k).lt.xr(i,1)) goto 6
      if(xwc(i,j,k).gt.xr(i,2)) goto 6
    7 continue
      new=new+1
      do 8 l=1,3
      do 9 j=1,nwa
      xwc(l,j,new)=xwc(l,j,k)
      vwc(l,j,new)=vwc(l,j,k)
    9 continue
    8 continue
    6 continue
      nwmc=new
      endif
    1 continue
      else
      new=0
      do 10 j=1,nwm
      idel=0
      do 11 i=1,3
      if(xr(i,1).lt.xr(i,2)) then
      do 12 k=1,nwa
      if(xw(i,k,j).gt.xr(i,1).and.xw(i,k,j).lt.xr(i,2)) idel=idel+1
   12 continue
      else
      idel=idel+nwa
      endif
   11 continue
      if(idel.ne.3*nwa) then
      new=new+1
      do 13 l=1,3
      do 14 k=1,nwa
      xw(l,k,new)=xw(l,k,j)
      vw(l,k,new)=vw(l,k,j)
   14 continue
   13 continue
      endif
   10 continue
      nwm=new
      new=0
      do 15 j=1,nwmc
      idel=0
      do 16 i=1,3
      if(xr(i,1).lt.xr(i,2)) then
      do 17 k=1,nwa
      if(xwc(i,k,j).gt.xr(i,1).and.xwc(i,k,j).lt.xr(i,2)) idel=idel+1
   17 continue
      else
      idel=idel+nwa
      endif
   16 continue
      if(idel.ne.3*nwa) then
      new=new+1
      do 18 l=1,3
      do 19 k=1,nwa
      xwc(l,k,new)=xwc(l,k,j)
      vwc(l,k,new)=vwc(l,k,j)
   19 continue
   18 continue
      endif
   15 continue
      nwmc=new
      endif
c
      pre_removw=.true.
c
      return
      end
      
