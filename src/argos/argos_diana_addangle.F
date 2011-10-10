      subroutine argos_diana_addangle()
c
c     add angle to the angle list
c     -------------------------------
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      logical lequiv
      integer last
      integer i,idangl,nanglp
      integer ndims(4)
c
      read(card(8:13),1000) lequiv,idangl
 1000 format(l1,i5)
c
      nanglp=nangl
      nangl=nangl+1
      if(nangl.gt.mxangl) call md_abort('Increase dimension mxangl',0)
c
      if(me.eq.0) then
      do 1 i=1,3
      read(lfncmd,1001,end=9999) card
 1001 format(a)
      call argos_diana_index(byte_mb(i_snam),
     + iangl(nangl,i,1),iangl(nangl,i,2))
    1 continue
      iangl(nangl,1,3)=idangl
      iangl(nangl,2,3)=1
c
      last=nangl
c      if(lequiv) then
c      call argos_diana_equitors(int_mb(i_iequi),int_mb(i_imol),
c     + idtors)
c      endif
      call argos_diana_listangl(byte_mb(i_snam),last)
c
      iangl(nangl,2,3)=nangl-nanglp
      endif
c
      call ga_brdcst(mag_d09,iangl,
     + 12*mxangl*ma_sizeof(mt_int,1,mt_byte),0)
c
      ndims(1)=0
      ndims(2)=0
      ndims(3)=nangl
      ndims(4)=0
c
      call ga_brdcst(mag_d10,ndims,4*ma_sizeof(mt_int,1,mt_byte),0)
c
      nangl=ndims(3)
c
      return
 9999 continue
      call md_abort('Failed to read angle',0)
      return
      end
c $Id$
