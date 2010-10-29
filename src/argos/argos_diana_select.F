      subroutine argos_diana_select(sgmnam,isel,wt,x,numsel)
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
      integer numsel
c
      integer i,j
      real*8 rtag,rtag2,r2
      character*80 target
      integer l_list,i_list
c
c     do selection on node 0
c
      if(me.eq.0) then
c
      read(card(11:80),1000) rtag,target
 1000 format(f12.6,a)
c
      if(.not.ma_push_get(mt_int,nsa,'list',l_list,i_list))
     + call md_abort('Failed to allocate list',0)
c
      call argos_diana_selection(sgmnam,nsa,target,int_mb(i_list))
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
