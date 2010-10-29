      logical function pre_removw(ir,xr,
     + xw,vw,mwm,mwa,nwm,nwa,xwc,vwc,mwmc,nwmc,xs,msa,nsa)
c
c $Id$
c
      implicit none
c
      integer ir
      real*8 xr(3,2)
      integer mwa,mwm,nwa,nwm,mwmc,nwmc,msa,nsa
      real*8 xw(3,mwa,mwm),xwc(3,mwa,mwmc)
      real*8 vw(3,mwa,mwm),vwc(3,mwa,mwmc)
      real*8 xs(3,msa)
c
      integer i,j,k,l,idel,new
      real*8 x2,d2
c
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
      elseif(ir.eq.1) then
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
      elseif(ir.eq.2) then
      new=0
      x2=xr(1,1)**2
      do 20 j=1,nwm
      idel=1
      do 21 k=1,nwa
      do 22 i=1,nsa
      d2=(xs(1,i)-xw(1,k,j))**2+(xs(2,i)-xw(2,k,j))**2+
     + (xs(3,i)-xw(3,k,j))**2
      if(d2.lt.x2) then
      idel=0
      goto 23
      endif
   22 continue
   21 continue
   23 continue
      if(idel.eq.0) then
      new=new+1
      do 24 k=1,nwa
      xw(1,k,new)=xw(1,k,j)
      xw(2,k,new)=xw(2,k,j)
      xw(3,k,new)=xw(3,k,j)
      vw(1,k,new)=vw(1,k,j)
      vw(2,k,new)=vw(2,k,j)
      vw(3,k,new)=vw(3,k,j)
   24 continue
      endif
   20 continue
      nwm=new
      new=0
      do 25 j=1,nwmc
      idel=1
      do 26 k=1,nwa
      do 27 i=1,nsa
      d2=(xs(1,i)-xwc(1,k,j))**2+(xs(2,i)-xwc(2,k,j))**2+
     + (xs(3,i)-xwc(3,k,j))**2
      if(d2.lt.x2) then
      idel=0
      goto 28
      endif
   27 continue
   26 continue
   28 continue
      if(idel.eq.0) then
      new=new+1
      do 29 k=1,nwa
      xwc(1,k,new)=xwc(1,k,j)
      xwc(2,k,new)=xwc(2,k,j)
      xwc(3,k,new)=xwc(3,k,j)
      vwc(1,k,new)=vwc(1,k,j)
      vwc(2,k,new)=vwc(2,k,j)
      vwc(3,k,new)=vwc(3,k,j)
   29 continue
      endif
   25 continue
      nwmc=new
      endif
c
      pre_removw=.true.
c
      return
      end
      
