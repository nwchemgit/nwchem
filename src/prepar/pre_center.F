      logical function pre_center(xw,mwm,mwa,nwm,nwa,xs,msa,nsa)
c
      implicit none
c
      integer mwa,mwm,nwa,nwm,msa,nsa
      real*8 xw(3,mwa,mwm),xs(3,msa)
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
      write(*,'(i5,f12.6)') k,s
      do 3 i=1,nsa
      xs(k,i)=xs(k,i)-s
    3 continue
      if(nwm.gt.0) then
      do 4 i=1,nwm
      do 5 j=1,nwa
      xw(k,j,i)=xw(k,j,i)-s
    5 continue
    4 continue
      endif
    1 continue
c
      pre_center=.true.
      return
      end
      
