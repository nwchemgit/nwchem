      subroutine argos_diana_define(sgmnam,isel,wt,x)
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*16 sgmnam(nsa)
      real*8 wt(nsa)
      integer isel(nsa)
      real*8 x(nsa,3)
c
      integer i,k,l
      integer id
      real*8 rtag,rtag2,r2
      character*80 target
      integer l_list,i_list
c
c     do selection on node 0
c
      if(me.eq.0) then
c
      read(card(8:80),1000) id,rtag,target
 1000 format(i7,f12.6,a)
c
      if(.not.ma_push_get(mt_int,nsa,'list',l_list,i_list))
     + call md_abort('Failed to allocate list',0)
c
      call argos_diana_selection(sgmnam,nsa,target,int_mb(i_list))
c
      if(abs(id).gt.mxdef) call md_abort('Group number too large',id)
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
      endif
c
      do 3 i=1,nsa
      if(int_mb(i_list-1+i).eq.1) then
      do 4 k=1,ldef(id)
      if(idef(id,k).eq.i) goto 3
    4 continue
      ldef(id)=ldef(id)+1
      if(ldef(id).gt.mxnum) call md_abort('Too many atoms in group',id)
      idef(id,ldef(id))=i
      if(rtag.gt.0.0d0) then
      rtag2=rtag*rtag
      do 5 k=1,nsa
      r2=(x(k,1)-x(i,1))**2+(x(k,2)-x(i,2))**2+(x(k,3)-x(i,3))**2
      if(r2.le.rtag2) then
      do 6 l=1,ldef(id)
      if(idef(id,l).eq.k) goto 5
    6 continue
      ldef(id)=ldef(id)+1
      if(ldef(id).gt.mxnum) call md_abort('Too many atoms in group',id)
      idef(id,ldef(id))=k
    7 continue
      endif
    5 continue
      endif
      endif
    3 continue
c
      if(.not.ma_pop_stack(l_list))
     + call md_abort('Failed to deallocate list',0)
c
      write(*,1001) id,ldef(id)
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
