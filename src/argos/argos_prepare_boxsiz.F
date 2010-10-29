      logical function argos_prepare_boxsiz(xw,xwc,mwm,mwmc,mwa,nwm,
     + nwmc,nwa,
     + xs,msa,nsa,box,iopt)
c
c $Id$
c
      implicit none
c
      integer mwm,mwmc,mwa,nwm,nwmc,nwa,msa,nsa,iopt
      real*8 xw(3,mwa,mwm),xwc(3,mwa,mwmc),xs(3,msa),box(3)
c
      integer i,j,k
      real*8 s
c
      if(nsa.le.0.and.nwm.le.0.and.nwmc.le.0) goto 9999
c
      do 1 k=1,3
      s=0.0d0
      do 2 i=1,nsa
      if(xs(k,i).gt.s) s=xs(k,i)
      if(xs(k,i).lt.-s) s=-xs(k,i)
    2 continue
      if(nwm.gt.0) then
      do 3 i=1,nwm
      do 4 j=1,nwa
      if(xw(k,j,i).gt.s) s=xw(k,j,i)
      if(xw(k,j,i).lt.-s) s=-xw(k,j,i)
    4 continue
    3 continue
      endif
      if(nwmc.gt.0) then
      do 5 i=1,nwmc
      do 6 j=1,nwa
      if(xwc(k,j,i).gt.s) s=xwc(k,j,i)
      if(xwc(k,j,i).lt.-s) s=-xwc(k,j,i)
    6 continue
    5 continue
      endif
      box(k)=s+s
    1 continue
c
      if(iopt.eq.1.or.iopt.eq.4) then
      s=0.0d0
      do 7 i=1,nsa
      s=max(s,xs(2,i)**2+xs(3,i)**2)
    7 continue
      if(nwm.gt.0) then
      do 8 i=1,nwm
      do 9 j=1,nwa
      s=max(s,xw(2,j,i)**2+xw(3,j,i)**2)
    9 continue
    8 continue
      endif
      if(nwmc.gt.0) then
      do 10 i=1,nwmc
      do 11 j=1,nwa
      s=max(s,xwc(2,j,i)**2+xwc(3,j,i)**2)
   11 continue
   10 continue
      endif
      box(2)=2.0d0*sqrt(s)
      box(3)=2.0d0*sqrt(s)
      endif
      if(iopt.eq.2.or.iopt.eq.4) then
      s=0.0d0
      do 12 i=1,nsa
      s=max(s,xs(1,i)**2+xs(3,i)**2)
   12 continue
      if(nwm.gt.0) then
      do 13 i=1,nwm
      do 14 j=1,nwa
      s=max(s,xw(1,j,i)**2+xw(3,j,i)**2)
   14 continue
   13 continue
      endif
      if(nwmc.gt.0) then
      do 15 i=1,nwmc
      do 16 j=1,nwa
      s=max(s,xwc(1,j,i)**2+xwc(3,j,i)**2)
   16 continue
   15 continue
      endif
      box(1)=2.0d0*sqrt(s)
      box(3)=2.0d0*sqrt(s)
      endif
      if(iopt.eq.3.or.iopt.eq.4) then
      s=0.0d0
      do 17 i=1,nsa
      s=max(s,xs(1,i)**2+xs(2,i)**2)
   17 continue
      if(nwm.gt.0) then
      do 18 i=1,nwm
      do 19 j=1,nwa
      s=max(s,xw(1,j,i)**2+xw(2,j,i)**2)
   19 continue
   18 continue
      endif
      if(nwmc.gt.0) then
      do 20 i=1,nwmc
      do 21 j=1,nwa
      s=max(s,xwc(1,j,i)**2+xwc(2,j,i)**2)
   21 continue
   20 continue
      endif
      box(1)=2.0d0*sqrt(s)
      box(2)=2.0d0*sqrt(s)
      endif
c
      argos_prepare_boxsiz=.true.
      return
c
 9999 continue
      box(1)=0.0d0
      box(2)=0.0d0
      box(3)=0.0d0
      argos_prepare_boxsiz=.true.
      return
      end
