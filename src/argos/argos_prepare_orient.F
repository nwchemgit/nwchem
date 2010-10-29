      logical function argos_prepare_orient(lfnout,xw,mwm,mwa,nwm,nwa,
     + xwc,mwmc,nwmc,xs,msa,nsa)
c
c $Id$
c
      implicit none
c
#include "util.fh"
c
      integer mwm,mwa,nwm,nwa,mwmc,nwmc,msa,nsa,lfnout
      real*8 xw(3,mwa,mwm),xwc(3,mwa,mwmc),xs(3,msa)
c
      integer i,ismax
      real*8 angle,dist,o(3),p(3),rsmax
c
      if(util_print('restart',print_low)) then
      write(lfnout,1000)
 1000 format(/,' Orientation solute')
      endif
c
      if(nsa.le.0) then
      argos_prepare_orient=.false.
      return
      endif
      if(nsa.le.1) then
      argos_prepare_orient=.true.
      return
      endif
c
c     find largest solute distance from origin
c
      ismax=0
      rsmax=0.0d0
      do 1 i=1,nsa
      dist=xs(1,i)**2+xs(2,i)**2+xs(3,i)**2
      if(dist.gt.rsmax) then
      ismax=i
      rsmax=dist
      endif
    1 continue
c
c     find rotation axis
c
      p(1)=xs(2,ismax)
      p(2)=-xs(1,ismax)
      p(3)=0.0d0
c
c     find rotation angle
c
      angle=acos(xs(3,ismax)/sqrt(rsmax))
c
c     find offset
c
      o(1)=0.0d0
      o(2)=0.0d0
      o(3)=0.0d0
c
c     rotate
c
      call argos_prepare_rotate(o,p,angle,xw,mwm,mwa,nwm,nwa,xwc,mwmc,
     + nwmc,
     + xs,msa,nsa)
      if(nsa.le.2) then
      argos_prepare_orient=.true.
      return
      endif
c
c     find largest distance in xy plane
c
      ismax=0
      rsmax=0.0d0
      do 2 i=1,nsa
      dist=xs(1,i)**2+xs(2,i)**2
      if(dist.gt.rsmax) then
      ismax=i
      rsmax=dist
      endif
    2 continue
c
c     find rotation axis
c
      p(1)=0.0d0
      p(2)=0.0d0
      p(3)=1.0d0
c
c     find rotation angle
c
      angle=acos(xs(2,ismax)/sqrt(rsmax))
c
c     find offset
c
      o(1)=0.0d0
      o(2)=0.0d0
      o(3)=0.0d0
c
c     rotate
c
      call argos_prepare_rotate(o,p,angle,xw,mwm,mwa,nwm,nwa,xwc,mwmc,
     + nwmc,
     + xs,msa,nsa)
c
      argos_prepare_orient=.true.
      return
      end
