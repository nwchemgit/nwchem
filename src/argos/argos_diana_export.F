      subroutine argos_diana_export()
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "argos_diana_common.fh"
#include "mafdecls.fh"
c
      integer i,j,ilt,iht,jlt,jht,nfr
      integer id1,id2,idexp,jtors
      logical lopen
      integer ntdim
c
      ntdim=ntval+1
      if(ntval.le.0) return
c
      call ga_sync()
c
      if(me.eq.0) then
c
      inquire(unit=lfndat,opened=lopen)
      if(.not.lopen)
     + open(unit=lfndat, file=fildat(1:index(fildat,' ')-1),
     + form='formatted',status='unknown')
c
      write(lfndat,1000) ntdim,nfrtot
 1000 format(i5,i10,/,'Time')
      if(ndist.gt.0) then
      do 1 i=1,ndist
      write(lfndat,1001) 'Distance ',cdist(i)
 1001 format(a,a)
    1 continue
      endif
      if(nangl.gt.0) then
      do 2 i=1,nangl
      write(lfndat,1002) 'Angle    ',cangl(i)
 1002 format(a,a)
    2 continue
      endif
      if(ntors.gt.0) then
      do 3 i=1,ntors
      write(lfndat,1003) 'Torsion  ',ctors(i)
 1003 format(a,a)
    3 continue
      endif
      if(nimps.gt.0) then
      endif
c
      if(.not.ma_push_get(mt_dbl,nfrdim*ntdim,'tv',l_tval,i_tval))
     + call md_abort('Failed to allocate tv',0)
c
      do 6 i=0,np-1
      call ga_distribution(ga_tval,i,ilt,iht,jlt,jht)
      call ga_get(ga_tval,ilt,iht,jlt,jht,dbl_mb(i_tval+1),ntdim)
      call ga_distribution(ga_trt,i,ilt,iht,jlt,jht)
      call ga_get(ga_trt,1,1,jlt,jht,dbl_mb(i_tval),ntdim)
      nfr=min(nfrdim,nfrtot-i*nfrdim)
      call argos_diana_wrtdat(dbl_mb(i_tval),ntdim,nfr)
    6 continue
c
      if(.not.ma_pop_stack(l_tval))
     + call md_abort('Failed to deallocate tv',0)
c
      endif
c
      return
      end
c $Id$
