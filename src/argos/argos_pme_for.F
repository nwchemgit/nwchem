      subroutine argos_pme_for(fw,nwmloc,fs,isq,isgm,nsaloc,
     + lnode,kfrom,kto)
c
      implicit none
c
#include "argos_pme_common.fh"
#include "mafdecls.fh"
c
      integer nwmloc,nsaloc
      real*8 fw(mwm,3,mwa),fs(msa,3)
      integer isq(msa),isgm(msa)
      integer lnode(np),kfrom(np),kto(np)
c
      integer inode,node,ngznon
c
c     wait until potential grid is available
c
      call timer_start(39)
      call argos_pme_wait(0,1)
      call timer_stop(39)
      call timer_start(40)
c
      do 1 inode=1,np
      node=me+inode
      if(node.gt.np) node=node-np
      node=node-1
c
c     node can be skipped if it did not contribute to the charge grid
c
      if(lnode(node+1).eq.0) goto 1
c
c     determine the number of planes on non-local node
c
      ngznon=0
      if(kto(node+1).gt.0) ngznon=kto(node+1)-kfrom(node+1)+1
      if(ngznon.eq.0) goto 1
c
      call timer_stop(40)
      call timer_start(41)
      call argos_pme_getgrd(node,ngznon,log_mb(i_lsld),int_mb(i_slmax),
     + dbl_mb(i_grd))
      call timer_stop(41)
      call timer_start(40)
c
      call argos_cafe_fpme(fw,nwmloc,fs,isq,isgm,nsaloc,
     + dbl_mb(i_grd),mgz,
     + int_mb(i_knod),int_mb(i_kfr),log_mb(i_lsld),node)
c
    1 continue
c
      call timer_stop(40)
c
      return
      end
c $Id$
