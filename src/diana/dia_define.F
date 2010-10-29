      subroutine dia_define(sgmnam,isel,wt,x)
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
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
      call dia_selection(sgmnam,nsa,target,int_mb(i_list))
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
      subroutine dia_wdefin()
c
      implicit none
c
#include "dia_common.fh"
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
      subroutine dia_selection(atomname,natoms,target,list)
c
      implicit none
c
      integer natoms
      character*16 atomname(natoms)
      character*80 target
      character*10 string,aname
      integer list(natoms)
c
      integer i,j,ip,ndxcol,length,isgmfr,isgmto,isegment
      logical lsgm,latm
c
      do 1 i=1,natoms
      list(i)=0
 1    continue
c
      length=index(target,' ')-1
      ndxcol=index(target,':')
      lsgm=ndxcol.gt.0
      latm=ndxcol.gt.0
      if(ndxcol.eq.0) then
      latm=.false.
      do 2 i=1,length
      lsgm=.false.
      if(target(i:i).eq.'1') lsgm=.true.
      if(target(i:i).eq.'2') lsgm=.true.
      if(target(i:i).eq.'3') lsgm=.true.
      if(target(i:i).eq.'4') lsgm=.true.
      if(target(i:i).eq.'5') lsgm=.true.
      if(target(i:i).eq.'6') lsgm=.true.
      if(target(i:i).eq.'7') lsgm=.true.
      if(target(i:i).eq.'8') lsgm=.true.
      if(target(i:i).eq.'9') lsgm=.true.
      if(target(i:i).eq.'0') lsgm=.true.
      if(target(i:i).eq.',') lsgm=.true.
      if(target(i:i).eq.'-') lsgm=.true.
      if(.not.lsgm) latm=.true.
      if(latm) goto 3
 2    continue
 3    continue
      if(lsgm) then
      ndxcol=length+1
      else
      ndxcol=0
      endif
      endif
c
      if(lsgm) then
      string='          '
      isgmto=0
      isgmfr=0
      ip=1
      do 4 i=1,ndxcol
      if(target(i:i).eq.','.or.i.eq.ndxcol) then
      string(10+ip-i:10)=target(ip:i-1)
      read(string,'(i10)') isgmto
      string='          '
      ip=i+1
      if(isgmfr.eq.0) isgmfr=isgmto
      elseif(target(i:i).eq.'-') then
      string(10+ip-i:10)=target(ip:i-1)
      read(string,'(i10)') isgmfr
      string='          '
      ip=i+1
      endif
      if(isgmto.gt.0) then
      do 5 j=1,natoms
      read(atomname(j)(12:16),'(i5)') isegment
      if(isegment.ge.isgmfr.and.isegment.le.isgmto) list(j)=-1
 5    continue
      isgmto=0
      isgmfr=0
      endif
 4    continue
      else
      if(latm) then
      do 6 j=1,natoms
      list(j)=-1
 6    continue
      else
      do 7 j=1,natoms
      list(j)=1
 7    continue
      endif
      endif
c
      if(latm) then
      string='          '
      ip=1
      do 8 i=ndxcol+1,length+1
      if(target(i:i).eq.','.or.i.eq.length+1) then
      do 9 j=1,natoms
      if(list(j).eq.-1) then
      aname(1:6)=atomname(j)(6:11)
      if(string(1:1).eq.'_') aname(1:1)='_'
      if(string(1:1).eq.'?') aname(1:1)='?'
      if(string(2:2).eq.'?') aname(2:2)='?'
      if(string(3:3).eq.'?') aname(3:3)='?'
      if(string(4:4).eq.'?') aname(4:4)='?'
      if(string(5:5).eq.'?') aname(5:5)='?'
      if(string(6:6).eq.'?') aname(6:6)='?'
      if(string(1:6).eq.aname(1:6)) list(j)=1
      endif
 9    continue
      string='          '
      ip=1
      else
      string(ip:ip)=target(i:i)
      ip=ip+1
      endif
 8    continue
      else
      do 10 i=1,natoms
      if(list(i).lt.0) list(i)=1
   10 continue
      endif
c
      do 11 i=1,natoms
      if(list(i).lt.0) list(i)=0
   11 continue
c
      return
      end
