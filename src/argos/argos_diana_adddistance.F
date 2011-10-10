      subroutine argos_diana_adddistance()
c
c     add distance to the distance list
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
      integer i,iddist,ndistp
      integer ndims(4)
c
      read(card(8:13),1000) lequiv,iddist
 1000 format(l1,i5)
c
      ndistp=ndist
      ndist=ndist+1
      if(ndist.gt.mxdist) call md_abort('Increase dimension mxdist',0)
c
      if(me.eq.0) then
      do 1 i=1,2
      read(lfncmd,1001,end=9999) card
 1001 format(a)
      call argos_diana_index(byte_mb(i_snam),
     + idist(ndist,i,1),idist(ndist,i,2))
    1 continue
      idist(ndist,1,3)=iddist
      idist(ndist,2,3)=1
c
      last=ndist
c      if(lequiv) then
c      call argos_diana_equitors(int_mb(i_iequi),int_mb(i_imol),
c     + idtors)
c      endif
      call argos_diana_listdist(byte_mb(i_snam),last)
c
      idist(ndist,2,3)=ndist-ndistp
      endif
c
      call ga_brdcst(mag_d09,idist,
     + 12*mxdist*ma_sizeof(mt_int,1,mt_byte),0)
c
      ndims(1)=0
      ndims(2)=0
      ndims(3)=ndist
      ndims(4)=0
c
      call ga_brdcst(mag_d10,ndims,4*ma_sizeof(mt_int,1,mt_byte),0)
c
      ndist=ndims(3)
c
      return
 9999 continue
      call md_abort('Failed to read distance',0)
      return
      end
c $Id$
