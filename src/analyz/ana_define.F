      subroutine ana_define(card,sgmnam,isel,wt,x)
c
c $Id: ana_define.F,v 1.1 2000-08-31 16:07:42 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
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
      integer i,j,k,l
      integer id,is,js,it
      character*5 cs,ct
      real*8 rtag,rtag2,r2
      character*80 card
c
c     do selection on node 0
c
      if(me.eq.0) then
c
      read(card(8:45),1000) id,is,js,rtag,cs
 1000 format(3i7,f12.6,a5)
c
      if(abs(id).gt.mxdef) call errquit('Group number too large',id)
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
      if(id.gt.ndef) call errquit('Groups to be defined in order',id)
c
      if(cs(1:1).eq.'_') cs(1:1)=' '
      do 3 i=1,nsa
      read(sgmnam(i),'(5x,a5,i6)') ct,it
      do 4 j=1,5
      if(cs(j:j).eq.'?') ct(j:j)='?'
    4 continue
      if((is.eq.0.or.(is.le.it.and.js.ge.it)).and.
     + (cs.eq.'     '.or.cs.eq.ct)) then
      do 6 k=1,ldef(id)
      if(idef(id,k).eq.i) goto 7
    6 continue
      ldef(id)=ldef(id)+1
      if(ldef(id).gt.mxnum) call errquit('Too many atoms in group',id)
      idef(id,ldef(id))=i
    7 continue
      if(rtag.gt.0.0d0) then
      rtag2=rtag*rtag
      do 5 k=1,nsa
      r2=(x(k,1)-x(i,1))**2+(x(k,2)-x(i,2))**2+(x(k,3)-x(i,3))**2
      if(r2.le.rtag2) then
      do 8 l=1,ldef(id)
      if(idef(id,l).eq.i) goto 9
    8 continue
      ldef(id)=ldef(id)+1
      if(ldef(id).gt.mxnum) call errquit('Too many atoms in group',id)
      idef(id,ldef(id))=i
    9 continue
      endif
    5 continue
      endif
      endif
    3 continue
      endif
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
