      subroutine ana_select(card,sgmnam,isel,wt,x,numsel)
c
c $Id$
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
      integer numsel
c
      integer i,j
      real*8 rtag,rtag2,r2
      character*80 card,target
      integer l_list,i_list
c
c     do selection on node 0
c
      if(me.eq.0) then
c
      read(card(8:80),1000) rtag,target
 1000 format(f12.6,a)
c
      if(.not.ma_push_get(mt_int,nsa,'list',l_list,i_list))
     + call md_abort('Failed to allocate list',0)
c
      call ana_selection(sgmnam,nsa,target,int_mb(i_list))
c
      do 1 i=1,nsa
      if(int_mb(i_list-1+i).eq.1) then
      isel(i)=1
      if(rtag.gt.0.0d0) then
      rtag2=rtag*rtag
      do 2 j=1,nsa
      r2=(x(j,1)-x(i,1))**2+(x(j,2)-x(i,2))**2+(x(j,3)-x(i,3))**2
      if(r2.le.rtag2) isel(j)=1
    2 continue
      endif
      endif
    1 continue
c
      if(.not.ma_pop_stack(l_list))
     + call md_abort('Failed to deallocate list',0)
c
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d06,isel,nsa*ma_sizeof(mt_int,1,mt_byte),0)
      endif
c
      numsel=0
      do 6 i=1,nsa
      wt(i)=1.0d0
      if(isel(i).eq.0) wt(i)=0.0d0
      numsel=numsel+isel(i)
    6 continue
c
      return
      end
      subroutine ana_all(isel,ival,wt,numsel)
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
      integer i,ival,numsel
c
      do 1 i=1,nsa
      isel(i)=ival
    1 continue
c
      numsel=0
      if(ival.ne.0) numsel=nsa
      do 2 i=1,nsa
      wt(i)=1.0d0
      if(isel(i).eq.0) wt(i)=0.0d0
    2 continue
c
      return
      end
      subroutine ana_wall(iwsel,ival,ndxw)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer iwsel(mwm),ndxw(mwm)
      integer i,ival
c
      nwsel=0
      if(ival.ne.0) nwsel=nwm
      do 1 i=1,nwm
      iwsel(i)=ival
      ndxw(i)=i
    1 continue
c
      return
      end
      subroutine ana_wselect(isel,xs,iwsel,xw,range,ndxw)
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer isel(msa),iwsel(mwm),ndxw(mwm)
      real*8 xs(msa,3),xw(mwm,mwa,3)
c
      integer i,j,k
      real*8 range,range2
c
      range2=range*range
c
c      if(me.eq.0) then
      do 1 i=1,nwm
      if(iwsel(i).eq.0) then
      do 2 j=ifr,ito
      if(isel(j).ne.0) then
      do 3 k=1,nwa
      if((xw(i,k,1)-xs(j,1))**2+(xw(i,k,2)-xs(j,2))**2+
     + (xw(i,k,3)-xs(j,3))**2.lt.range2) then
      iwsel(i)=1
      goto 1
      endif
    3 continue
      endif
    2 continue
      endif
    1 continue
      nwsel=0
      do 4 i=1,nwm
      if(iwsel(i).gt.0) then
      nwsel=nwsel+1
      ndxw(nwsel)=i
      endif
    4 continue
c
      if(np.gt.1) then
      call ga_igop(mag_d14,nwsel,1,'+')
      endif
c      endif
c
      return
      end
      subroutine ana_wcount(iwsel,ndxw)
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer iwsel(mwm),ndxw(mwm),i
c
      nwsel=0
      do 1 i=1,nwm
      if(iwsel(i).gt.0) then
      nwsel=nwsel+1
      ndxw(nwsel)=i
      endif
    1 continue
c
      return
      end
      subroutine ana_savsel(lfn,fil,iwsel,isel)
c
      implicit none
c
#include "ana_common.fh"
c
      integer lfn,iwsel(mwm),isel(msa)
      character*255 fil
c
      integer i
c
      if(me.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + form='formatted',status='unknown')
      rewind(lfn)
      write(lfn,1000) (iwsel(i),i=1,nwm)
      write(lfn,1000) (isel(i),i=1,nsa)
 1000 format(40i2)
      close(unit=lfn)
      write(*,'(/,a,a)') ' Selection saved to ',fil(1:index(fil,' ')-1)
      endif
c
      return
      end
      subroutine ana_readsel(lfn,fil,iwsel,isel)
c
      implicit none
c
#include "ana_common.fh"
c
      integer lfn,iwsel(mwm),isel(msa)
      character*255 fil
c
      integer i
c
      if(me.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + form='formatted',status='old')
      rewind(lfn)
      read(lfn,1000) (iwsel(i),i=1,nwm)
      read(lfn,1000) (isel(i),i=1,nsa)
 1000 format(40i2)
      close(unit=lfn)
      write(*,'(/,a,a)') ' Selection read from ',fil(1:index(fil,' ')-1)
      endif
c
      return
      end
