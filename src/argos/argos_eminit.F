      subroutine argos_eminit(xw,yw,xs,ys)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
c
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3)
c
      integer i,j
      real*8 dxmax
c
      if(nwmloc.gt.0) then
      do 1 j=1,nwa
      do 2 i=1,nwmloc
      yw(i,1,j)=xw(i,1,j)
      yw(i,2,j)=xw(i,2,j)
      yw(i,3,j)=xw(i,3,j)
    2 continue
    1 continue
      endif
c
      if(nsaloc.gt.0) then
      do 3 i=1,nsaloc
      ys(i,1)=xs(i,1)
      ys(i,2)=xs(i,2)
      ys(i,3)=xs(i,3)
    3 continue
      endif
c
      return
c
      lpair=.true.
      lload=.true.
      lhop=.false.
c
      call argos_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call argos_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
c     shake
c
      call argos_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
      return
      end
c $Id$
