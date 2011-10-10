      subroutine argos_diana_addtorsion()
c
c     add torsion to the torsion list
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
      integer i,idtors,ntorsp
      integer ndims(4)
c
      read(card(8:13),1000) lequiv,idtors
 1000 format(l1,i5)
c
      ntorsp=ntors
      ntors=ntors+1
      if(ntors.gt.mxtors) call md_abort('Increase dimension mxtors',0)
c
      if(me.eq.0) then
      do 1 i=1,4
      read(lfncmd,1001,end=9999) card
 1001 format(a)
      call argos_diana_index(byte_mb(i_snam),
     + itors(ntors,i,1),itors(ntors,i,2))
    1 continue
      itors(ntors,1,3)=idtors
      itors(ntors,2,3)=1
c
      last=ntors
      if(lequiv) then
      call argos_diana_equitors(int_mb(i_iequi),int_mb(i_imol),
     + idtors)
      endif
      call argos_diana_listtors(byte_mb(i_snam),last)
c
      itors(ntors,2,3)=ntors-ntorsp
      endif
c
      call ga_brdcst(mag_d09,itors,
     + 12*mxtors*ma_sizeof(mt_int,1,mt_byte),0)
c
      ndims(1)=0
      ndims(2)=0
      ndims(3)=ntors
      ndims(4)=0
c
      call ga_brdcst(mag_d10,ndims,4*ma_sizeof(mt_int,1,mt_byte),0)
c
      ntors=ndims(3)
c
      return
 9999 continue
      call md_abort('Failed to read torsion',0)
      return
      end
c $Id$
