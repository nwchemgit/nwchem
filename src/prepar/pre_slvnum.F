      logical function pre_slvnum(islv,boxs,nwms,box,nwm)
c
c $Id$
c
      implicit none
c
      integer islv,nwms,nwm
      real*8 boxs(3),box(3)
c
      nwm=0
      if(islv.eq.1) then
      nwm=2*int(dble(nwms)*(box(1)*box(2)*box(3))/
     + (boxs(1)*boxs(2)*boxs(3)))
      elseif(islv.eq.2) then
      nwm=int(dble(nwms)*0.5d0*(box(1)*box(2)*box(3))/
     + (boxs(1)*boxs(2)*boxs(3)))
      elseif(islv.eq.3) then
      nwm=2*int(dble(nwms)*(box(1)*box(2)*box(3))/
     + (boxs(1)*boxs(2)*boxs(3)))
      endif
c
      pre_slvnum=.true.
      return
      end
