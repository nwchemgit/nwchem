      logical function pre_boxsiz(xw,mwm,mwa,nwm,nwa,xs,msa,nsa,box)
c
      implicit none
c
      integer mwm,mwa,nwm,nwa,msa,nsa
      real*8 xw(3,mwa,mwm),xs(3,msa),box(3)
c
      integer i,j,k
      real*8 s
c
      if(nsa.le.0.and.nwm.le.0) goto 9999
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
      box(k)=s+s
    1 continue
c
      pre_boxsiz=.true.
      return
c
 9999 continue
      pre_boxsiz=.false.
      return
      end
