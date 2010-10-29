      subroutine argos_diana_center_ref(x)
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      real*8 x(nsa,3)
      integer i,j
c
c     center and fold
c
      if(me.eq.0.and.icent.gt.0) then
      if(.not.ltop) call md_abort('No topology for center',0)
      call argos_diana_center_fold(x,x,int_mb(i_imol),0)
      endif
c
c     broadcast coordinates to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d02,x(1,1),3*nsa*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
      do 1 j=1,3
      xsmin(j)=x(1,j)
      xsmax(j)=x(1,j)
      do 1 i=2,nsa
      xsmin(j)=min(xsmin(j),x(i,j))
      xsmax(j)=max(xsmax(j),x(i,j))
    1 continue
c
c      xsmin=xsmin-0.5d0
c      xsmax=xsmax+0.5d0
c
      return
      end
