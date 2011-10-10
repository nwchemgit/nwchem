      subroutine argos_diana_wdefin()
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer id
      character*5 cs,ct
      real*8 rtag
      integer i,j,k
c
      read(card(8:31),1000) id,rtag,cs
 1000 format(i7,f12.6,a5)
c
      if(id.eq.0) then
      ndef=0
      else
      if(id.lt.0) then
      id=-id
      ldef(id)=0
      endif
      if(id.eq.ndef+1) then
      ndef=id
      ldef(id)=0
      endif
      if(id.gt.ndef) call md_abort('Groups to be defined in order',id)
      if(ldef(id).gt.0) ldef(id)=0
c
      if(cs(1:1).eq.'_') cs(1:1)=' '
c
      do 3 i=1,nwa
      ct=wnam(i)(6:10)
      do 4 j=1,5
      if(cs(j:j).eq.'?') ct(j:j)='?'
    4 continue
      if(cs.eq.'     '.or.cs.eq.ct) then
      do 6 k=1,-ldef(id)
      if(idef(id,k).eq.i) goto 7
    6 continue
      ldef(id)=ldef(id)-1
      if(-ldef(id).gt.mxnum) call md_abort('Too many atoms in group',id)
      idef(id,-ldef(id))=i
    7 continue
      endif
    3 continue
c
      write(*,1001) id,-ldef(id)
 1001 format(' Atom group',i5,' contains ',i7,' atoms')
c
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d06,ndef,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mag_d06,ldef,ndef*ma_sizeof(mt_int,1,mt_byte),0)
      do 10 i=1,ndef
      call ga_brdcst(mag_d06,idef(i,1),
     + ldef(i)*ma_sizeof(mt_int,1,mt_byte),0)
   10 continue
      endif
c
      return
      end
c $Id$
