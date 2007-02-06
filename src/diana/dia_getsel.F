      subroutine dia_getsel
c
c $Id: dia_getsel.F,v 1.4 2007-02-06 21:24:16 d3j191 Exp $
c
c     get atom selection
c     ------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      logical lssup,lsout,lequiv
c
      read(card(8:10),1000) lssup,lsout,lequiv
 1000 format(3l1)
c
      if(lssup) then
      if(.not.lsels) then
      call dia_all(int_mb(i_ssel),0,dbl_mb(i_swt),nsel)
      lsels=.true.
      endif
      call dia_select(card,byte_mb(i_snam),int_mb(i_ssel),dbl_mb(i_swt),
     + dbl_mb(i_xref),nsel)
      if(me.eq.0) then
      write(*,1001) nsel,nsa
 1001 format(' Selected',i5,' out of',i5,' atoms for superposition')
      endif
      elseif(lsout) then
      if(.not.lselo) then
      call dia_all(int_mb(i_osel),0,dbl_mb(i_owt),nselo)
      lselo=.true.
      endif
      call dia_select(card,byte_mb(i_snam),int_mb(i_osel),dbl_mb(i_owt),
     + dbl_mb(i_xref),nselo)
      if(me.eq.0) then
      write(*,1002) nselo,nsa
 1002 format(' Selected',i5,' out of',i5,' atoms for output')
      endif
      else
      if(.not.lsel) then
      call dia_all(int_mb(i_isel),0,dbl_mb(i_wt),nsel)
      lsel=.true.
      endif
      call dia_select(byte_mb(i_snam),int_mb(i_isel),dbl_mb(i_wt),
     + dbl_mb(i_xref),nsel)
      if(me.eq.0) then
      write(*,1003) nsel,nsa
 1003 format(' Selected',i5,' out of',i5,' atoms')
      endif
      endif
c
      return
      end
      subroutine dia_getsls
c
c     get atom selection
c     ------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      if(.not.lsels) then
      call dia_all(int_mb(i_ssel),0,dbl_mb(i_swt),nsel)
      lsels=.true.
      endif
      call dia_select(card,byte_mb(i_snam),int_mb(i_ssel),dbl_mb(i_swt),
     + dbl_mb(i_xref),nsel)
      if(me.eq.0) then
      write(*,3533) nsel,nsa
 3533 format(' Selected',i5,' out of',i5,' atoms for superposition')
      endif
c
      return
      end
      subroutine dia_getslo
c
c     get atom selection
c     ------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      if(.not.lselo) then
      call dia_all(int_mb(i_osel),0,dbl_mb(i_owt),nselo)
      lselo=.true.
      endif
      call dia_select(card,byte_mb(i_snam),int_mb(i_osel),dbl_mb(i_owt),
     + dbl_mb(i_xref),nselo)
      if(me.eq.0) then
      write(*,3633) nselo,nsa
 3633 format(' Selected',i5,' out of',i5,' atoms for output')
      endif
c
      return
      end
      subroutine dia_addtorsion()
c
c     add torsion to the torsion list
c     -------------------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      logical lequiv
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
      call dia_index(byte_mb(i_snam),
     + itors(ntors,i,1),itors(ntors,i,2))
    1 continue
      itors(ntors,1,3)=idtors
      itors(ntors,2,3)=1
c
      if(lequiv) then
      call dia_equitors(byte_mb(i_snam),int_mb(i_iequi),int_mb(i_imol),
     + idtors)
      endif
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
      subroutine dia_equitors(sgmnam,iequi,imol,idtors)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      character*16 sgmnam(msa)
      integer iequi(msa),imol(msa)
      integer i,j,last,idtors
c
      last=ntors
c
      do 1 i=2,4
      if(imol(itors(ntors,i,1)).ne.imol(itors(ntors,1,1))) return
    1 continue
c
    2 continue
      do 3 i=1,4
      if(iequi(itors(ntors,i,1)).eq.0) goto 9
    3 continue
c
      ntors=ntors+1
      itors(ntors,1,3)=idtors
      if(ntors.gt.mxtors) call md_abort('Increase dimension mxtors',0)
c
      do 4 i=1,4
      itors(ntors,i,1)=iequi(itors(ntors-1,i,1))
      itors(ntors,i,2)=0
    4 continue
      goto 2
c
    9 continue
      if(me.eq.0) then
      do 5 i=last,ntors
      write(*,1000) i,(sgmnam(itors(i,j,1))(11:16),':',
     + sgmnam(itors(i,j,1))(1:10),j=1,4)
 1000 format(i4,4(a6,a,a10))
    5 continue
      endif
c
      return
      end
      subroutine dia_addorder()
c
c     add order parameter definition to the order parameter list
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      logical lequiv
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
      call dia_index(byte_mb(i_snam),
     + iord(nord,i,1),iord(nord,i,2))
    1 continue
      iord(nord,1,3)=idordr
      iord(nord,2,3)=1
      iord(nord,3,3)=ndordr
c
      if(lequiv) then
      call dia_equiord(byte_mb(i_snam),int_mb(i_iequi),int_mb(i_imol),
     + idordr,ndordr)
      endif
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
      subroutine dia_equiord(sgmnam,iequi,imol,idordr,ndordr)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      character*16 sgmnam(msa)
      integer iequi(msa),imol(msa)
      integer i,j,last,idordr,ndordr
c
      last=nord
c
      do 1 i=2,ndordr+1
      if(imol(iord(nord,i,1)).ne.imol(iord(nord,1,1))) return
    1 continue
c
    2 continue
      do 3 i=1,ndordr+1
      if(iequi(iord(nord,i,1)).eq.0) goto 9
    3 continue
c
      nord=nord+1
      iord(nord,1,3)=idordr
      iord(nord,3,3)=ndordr
      if(nord.gt.mord) call md_abort('Increase dimension mord',0)
c
      do 4 i=1,ndordr+1
      iord(nord,i,1)=iequi(iord(nord-1,i,1))
      iord(nord,i,2)=0
    4 continue
      goto 2
c
    9 continue
      if(me.eq.0) then
      do 5 i=last,nord
      write(*,1000) i,(sgmnam(iord(i,j,1))(11:16),':',
     + sgmnam(iord(i,j,1))(1:10),j=1,ndordr+1)
 1000 format(i4,4(a6,a,a10))
    5 continue
      endif
c
      return
      end
      
