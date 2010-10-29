      logical function pre_hnames(latm,catm,matm,natm,
     + lbnd,mbnd,nbnd)
c
c $Id$
c
      implicit none
c
      integer matm,natm,mbnd,nbnd
      integer latm(5,matm),lbnd(2,mbnd)
      character*6 catm(3,matm)
c
      integer i,j,nhb,ihb(10)
c
      do 1 i=1,natm
      nhb=0
      do 2 j=1,nbnd
      if(lbnd(1,j).eq.i.and.latm(2,lbnd(2,j)).eq.1) then
      nhb=nhb+1
      ihb(nhb)=lbnd(2,j)
      endif
      if(lbnd(2,j).eq.i.and.latm(2,lbnd(1,j)).eq.1) then
      nhb=nhb+1
      ihb(nhb)=lbnd(1,j)
      endif
    2 continue
      if(nhb.gt.0) then
      if(nhb.eq.1) then
      write(catm(2,ihb(1)),'(a2,a4)') ' H',catm(2,i)(3:6)
      else
      do 3 j=1,nhb
      write(catm(2,ihb(j)),'(i1,a1,a4)') j+1,'H',catm(2,i)(3:6)
    3 continue
      endif
      endif
    1 continue
c
      pre_hnames=.true.
c
      return
      end
