      subroutine qmd_forces(irtdb,isag,isat,xs,fs,msa,nsa,energy,lesp)
c
c $Id$
c
      implicit none
#include "errquit.fh"
c
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "qmd_common.fh"
#include "nwmd_constants.fh"
c
      external task_gradient,esp
      logical task_gradient,esp
c
c      data lgeom/0/
c
      integer irtdb
      integer msa,nsa
      integer isag(msa),isat(msa)
      real*8 xs(msa,3),fs(msa,3),energy
      logical lesp
c
      integer i,isa
      integer i_qmdq,l_qmdq,i_qmdt,l_qmdt,i_qmdb,l_qmdb
      real*8 c(3),q
      character*16 tag
c
      if(lgeom.eq.0) then
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmd: failed to create geometry',0, GEOM_ERR)
      if(.not.ma_push_get(mt_dbl,nsa,'qmdq',l_qmdq,i_qmdq))
     + call errquit('Failed to allocate qmdq',0, MA_ERR)
      if(.not.ma_push_get(mt_byte,16*nsa,'qmdt',l_qmdt,i_qmdt))
     + call errquit('Failed to allocate qmdt',0, MA_ERR)
      call qmd_geom_init(nsa,isag,isat,dbl_mb(i_qmdq),byte_mb(i_qmdt))
      if(ga_nodeid().eq.0) then
      if(.not.geom_cart_set(igeom,nsa,byte_mb(i_qmdt),xs,
     + dbl_mb(i_qmdq))) call errquit('Failed to initialize geom',0,
     &       GEOM_ERR)
      endif
      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit('Failed to store geom to rtdb',0, RTDB_ERR)
      if(.not.ma_pop_stack(l_qmdt))
     + call errquit('Failed to deallocate qmdt',0, MA_ERR)
      if(.not.ma_pop_stack(l_qmdq))
     + call errquit('Failed to deallocate qmdq',0, MA_ERR)
      lgeom=1
      endif
c
      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit('Failed to load geometry from rtdb',0, RTDB_ERR)
c
      do 1 i=1,nsa
      isa=isag(i)
      if(.not.geom_cent_get(igeom,isa,tag,c,q))
     + call errquit('geom_cent_get failed',isa, GEOM_ERR)
      c(1)=xs(isa,1)*cnm2au
      c(2)=xs(isa,2)*cnm2au
      c(3)=xs(isa,3)*cnm2au
      if(.not.geom_cent_set(igeom,i,tag,c,q))
     + call errquit('geom_cent_set failed',isa, GEOM_ERR)
    1 continue
      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit('geom_rtdb_store failed',0, RTDB_ERR)
c
      call ga_sync()
      if(jdebug.gt.0) call ma_summarize_allocated_blocks()
      call ga_sync()
      if(.not.task_gradient(irtdb))
     + call errquit('qmd: Task_gradient failed',0, CALC_ERR)
      call ga_sync()
      if(lesp) then
      if(.not.esp(irtdb)) call errquit('qmd esp error',0, RTDB_ERR)
      endif
c
      if(.not.ma_push_get(mt_dbl,3*msa,'qmdb',l_qmdb,i_qmdb))
     + call errquit('Failed to allocate qmdq',0, MA_ERR)
      if(.not.ma_verify_allocator_stuff())
     + call errquit('Error qmd_forces 1',0, MA_ERR)
      call qmd_f(irtdb,isag,dbl_mb(i_qmdb),fs,msa,nsa,energy)
      if(.not.ma_verify_allocator_stuff())
     + call errquit('Error qmd_forces 2',0, MA_ERR)
      if(.not.ma_pop_stack(l_qmdb))
     + call errquit('Failed to deallocate qmdb',0, MA_ERR)
c
      return
      end
      subroutine qmd_f(irtdb,isag,buf,fs,msa,nsa,energy)
c
      implicit none
#include "errquit.fh"
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "qmd_common.fh"
#include "nwmd_constants.fh"
c
      integer irtdb,msa,nsa
      integer isag(msa)
      real*8 buf(3,msa),fs(msa,3),energy
c
      integer i,isa
c
      if(.not.rtdb_get(irtdb,'task:gradient',mt_dbl,3*msa,buf))
     + call errquit('qmd_f: rtdb_get gradient failed',0, RTDB_ERR)
      if(.not.rtdb_get(irtdb,'task:energy',mt_dbl,1,energy))
     + call errquit('qmd_f: rtdb_get energy failed',0, RTDB_ERR)
c
      do 2 i=1,nsa
      isa=isag(i)
      fs(isa,1)=(-buf(1,i))*(cau2kj/cau2nm)
      fs(isa,2)=(-buf(2,i))*(cau2kj/cau2nm)
      fs(isa,3)=(-buf(3,i))*(cau2kj/cau2nm)
    2 continue
c
      energy=energy*cau2kj
c
      return
      end
      subroutine qmd_geom_init(nsa,isag,isat,q,tag)
c
      implicit none
c
      integer cf_element
      external cf_element
c
      integer nsa
      integer isag(nsa),isat(nsa)
      real*8 q(nsa)
      character*16 tag(nsa)
c
      integer i,ia
c
      do 1 i=1,nsa
      ia=cf_element(isat(i))
      call cf_num2tag(ia,tag(isag(i)))
      q(i)=dble(ia)
    1 continue
c
      return
      end
      subroutine qmd_wrttrj(lfntrj,mwm,nwm,mwa,nwa,msa,nsa,
     + lxw,lvw,lxs,lvs,box,stime,pres,temp,tempw,temps,xw,vw,xs,vs)
c
      integer lfntrj,msa,nsa,mwm,nwm,mwa,nwa
      logical lxw,lvw,lxs,lvs
      real*8 stime,pres,temp,tempw,temps,box(3)
      real*8 xw(mwm,3,mwa),vw(mwm,3,mwa)
      real*8 xs(msa,3),vs(msa,3)
c
      character*10 rdate,rtime
c
#include "global.fh"
c
      if(ga_nodeid().eq.0) then
c
      call swatch(rdate,rtime)
c
      write(lfntrj,1000)
 1000 format('frame')
      write(lfntrj,1001) stime,temp,pres,rdate,rtime
 1001 format(2f12.6,1pe12.5,1x,2a10)
      write(lfntrj,1002) box(1),0.0d0,0.0d0
      write(lfntrj,1002) 0.0d0,box(2),0.0d0
      write(lfntrj,1002) 0.0d0,0.0d0,box(3)
 1002 format(3f12.6)
      write(lfntrj,1003) lxw,lvw,lxs,lvs,nwm,nwa,nsa
 1003 format(4l1,3i10)
c
      if((lxw.or.lvw).and.nwm.gt.0) then
      if(lvw) then
      do 1 i=1,nwm
      write(lfntrj,1004) ((xw(i,j,k),j=1,3),(vw(i,j,k),j=1,3),k=1,nwa)
 1004 format(6f8.3)
    1 continue
      else
      do 2 i=1,nwm
      write(lfntrj,1005) ((xw(i,j,k),j=1,3),k=1,nwa)
 1005 format(3f8.3)
    2 continue
      endif
      endif
c
      if((lxs.or.lvs).and.nsa.gt.0) then
      if(lvs) then
      do 3 i=1,nsa
      write(lfntrj,1006) (xs(i,j),j=1,3),(vs(i,j),j=1,3)
 1006 format(6f8.3)
    3 continue
      else
      do 4 i=1,nsa
      write(lfntrj,1007) (xs(i,j),j=1,3)
 1007 format(3f8.3)
    4 continue
      endif
      endif
c
      endif
c
      return
      end
