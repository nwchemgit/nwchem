      subroutine ana_select(card,sgmnam,isel,wt,x)
c
c $Id: ana_select.F,v 1.5 2000-03-16 20:45:38 d3j191 Exp $
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
      integer i,j,k
      integer is,js,it
      character*5 cs,ct
      real*8 rtag,rtag2,r2
      character*80 card
c
c     do selection on node 0
c
      if(me.eq.0) then
c
      read(card(8:38),1000) is,js,rtag,cs
 1000 format(2i7,f12.6,a5)
      if(cs(1:1).eq.'_') cs(1:1)=' '
      do 3 i=1,nsa
      read(sgmnam(i),'(5x,a5,i6)') ct,it
      do 4 j=1,5
      if(cs(j:j).eq.'?') ct(j:j)='?'
    4 continue
      if((is.eq.0.or.(is.le.it.and.js.ge.it)).and.
     + (cs.eq.'     '.or.cs.eq.ct)) then
      isel(i)=1
      if(rtag.gt.0.0d0) then
      rtag2=rtag*rtag
      do 5 k=1,nsa
      r2=(x(k,1)-x(i,1))**2+(x(k,2)-x(i,2))**2+(x(k,3)-x(i,3))**2
      if(r2.le.rtag2) isel(k)=1
    5 continue
      endif
      endif
    3 continue
c
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d06,isel,nsa*ma_sizeof(mt_int,1,mt_byte),0)
      endif
c
      nsel=0
      do 6 i=1,nsa
      wt(i)=1.0d0
      if(isel(i).eq.0) wt(i)=0.0d0
      nsel=nsel+isel(i)
    6 continue
c
      return
      end
      subroutine ana_all(isel,ival,wt)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer isel(nsa)
      real*8 wt(nsa)
      integer i,ival
c
      do 1 i=1,nsa
      isel(i)=ival
    1 continue
c
      nsel=0
      do 2 i=1,nsa
      wt(i)=1.0d0
      if(isel(i).eq.0) wt(i)=0.0d0
      nsel=nsel+isel(i)
    2 continue
c
      return
      end
