      subroutine argos_diana_addorder()
c
c     add order parameter definition to the order parameter list
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      logical lequiv
      integer last
      integer i,idordr,ndordr,nordp
      integer ndims(4)
c
      read(card(8:18),1000) lequiv,idordr,ndordr
 1000 format(l1,2i5)
c
      nordp=nord
      nord=nord+1
      if(nord.gt.mord) call md_abort('Increase dimension mord',0)
c
      if(me.eq.0) then
      do 1 i=1,ndordr+1
      read(lfncmd,1001,end=9999) card
 1001 format(a)
      call argos_diana_index(byte_mb(i_snam),
     + iord(nord,i,1),iord(nord,i,2))
    1 continue
      iord(nord,1,3)=idordr
      iord(nord,2,3)=1
      iord(nord,3,3)=ndordr
c
      last=nord
      if(lequiv) then
      call argos_diana_equiord(int_mb(i_iequi),int_mb(i_imol),
     + idordr,ndordr)
      endif
      call argos_diana_listord(byte_mb(i_snam),last,ndordr)
c
      do 2 i=nordp+1,nord
      iord(i,2,3)=nord-nordp
    2 continue
      endif
c
      call ga_brdcst(mag_d09,iord,12*mord*ma_sizeof(mt_int,1,mt_byte),0)
c
      ndims(1)=0
      ndims(2)=0
      ndims(3)=nord
      ndims(4)=0
c
      call ga_brdcst(mag_d10,ndims,4*ma_sizeof(mt_int,1,mt_byte),0)
c
      nord=ndims(3)
c
      return
 9999 continue
      call md_abort('Failed to read order parameter',0)
      return
      end
c $Id$
