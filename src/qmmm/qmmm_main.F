      subroutine qmmm_forces(inrtdb,mwm,nwm,mwa,nwa,iwz,xw,fw,
     + msa,nsa,isat,isdt,isq,isz,xs,fs,energy)
c
c $Id: qmmm_main.F,v 1.1 2000-10-23 20:09:44 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "qmmm_common.fh"
c
      external task_gradient
      logical task_gradient
c
      integer inrtdb
      integer mwm,nwm,mwa,nwa,msa,nsa
      integer iwz(mwm),isat(msa),isdt(msa),isq(msa),isz(msa)
      real*8 energy,xw(mwm,3,mwa),fw(mwm,3,mwa),xs(msa,3),fs(msa,3)
c
      integer i,np,me,ioffq,ioffz,igeom
      integer i_ioff,l_ioff,i_c,l_c,i_q,l_q,i_t,l_t,i_i,l_i,i_m,l_m
      integer nzatom,nqatom,nzatomt,nqatomt
c
      np=ga_nnodes()
      me=ga_nodeid()
      irtdb=inrtdb
c
      if(.not.ma_push_get(mt_int,2*np,'ioff',l_ioff,i_ioff))
     + call errquit('qmmm: Failed to allocate memory for ioff',me)
c
      write(*,8888) me,(isdt(i),i=1,nsa)
 8888 format(/,' Solute dynamics types',/,1x,i1,':',25i3,/,(3x,25i3))
      write(*,8889) me,(isz(i),i=1,nsa)
 8889 format(/,' Solute zone types',/,1x,i1,'=',25i3,/,(3x,25i3))
c
      nzatom=0
      nqatom=0
      do 1 i=1,nsa
      if(iand(isdt(i),lquant).ne.0) nqatom=nqatom+1
      if(isz(i).gt.0) nzatom=nzatom+1
    1 continue
      do 2 i=1,nwm
      if(iwz(i).gt.0) nzatom=nzatom+mwa
    2 continue
c
      call qmmm_offset(me,np,nqatom,nzatom,ioffq,ioffz,nqatomt,nzatomt,
     + int_mb(i_ioff))
c
      if(.not.ma_pop_stack(l_ioff))
     + call errquit('qmmm: Failed to deallocate memory for ioff',me)
c
      if(.not.ma_push_get(mt_dbl,3*(nqatomt+nzatomt),'c',l_c,i_c))
     + call errquit('qmmm: Failed to allocate memory for c',me)
      if(.not.ma_push_get(mt_dbl,nqatomt+nzatomt,'q',l_q,i_q))
     + call errquit('qmmm: Failed to allocate memory for q',me)
      if(.not.ma_push_get(mt_dbl,nqatomt+nzatomt,'m',l_m,i_m))
     + call errquit('qmmm: Failed to allocate memory for m',me)
      if(.not.ma_push_get(mt_byte,16*(nqatomt+nzatomt),'t',l_t,i_t))
     + call errquit('qmmm: Failed to allocate memory for t',me)
      if(.not.ma_push_get(mt_int,nqatomt+nzatomt,'i',l_i,i_i))
     + call errquit('qmmm: Failed to allocate memory for i',me)
c
      call qmmm_geometry(nqatomt+nzatomt,dbl_mb(i_c),dbl_mb(i_q),
     + dbl_mb(i_m),byte_mb(i_t),int_mb(i_i),ioffq,ioffz,
     + mwm,nwm,mwa,msa,nsa,iwz,isat,isdt,isq,isz,xw,xs)
c
      if(.not.ma_pop_stack(l_i))
     + call errquit('qmmm: Failed to deallocate memory for i',me)
c
      if(.not.geom_create(igeom,'geometry'))
     + call errquit('qmmm: Failed to create geometry',0)
      if(.not.geom_cart_set(igeom,nqatomt+nzatomt,byte_mb(i_t),
     + dbl_mb(i_c),dbl_mb(i_q)))
     + call errquit('qmmm: Failed to initialize geometry',0)
      if(.not.geom_masses_set(igeom,nqatomt+nzatomt,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0)
      call geom_compute_values(igeom)
c
      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     + call errquit('qmmm: Failed to store geom to rtdb',0)
c
      if(.not.ma_pop_stack(l_t))
     + call errquit('qmmm: Failed to deallocate memory for t',me)
      if(.not.ma_pop_stack(l_m))
     + call errquit('qmmm: Failed to deallocate memory for m',me)
      if(.not.ma_pop_stack(l_q))
     + call errquit('qmmm: Failed to deallocate memory for q',me)
c
      call ga_sync()
      if(.not.task_gradient(irtdb))
     + call errquit('qmmm: Task_gradient failed',0)
      call ga_sync()
c
      call qmmm_f(nqatomt+nzatomt,dbl_mb(i_c),ioffq,ioffz,
     + mwm,nwm,mwa,msa,nsa,iwz,isat,isdt,isq,isz,fw,fs,energy)
c
      if(.not.ma_pop_stack(l_c))
     + call errquit('qmmm: Failed to deallocate memory for c',me)
c
      if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',me)
c
      return
      end
      subroutine qmmm_offset(me,np,nqa,nza,ioffq,ioffz,nqat,nzat,loff)
c
      implicit none
c
#include "qmmm_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer me,np,nqa,nza,nqat,nzat,ioffq,ioffz
      integer loff(np,2)
c
      integer i
c
      do 1 i=1,np
      loff(i,1)=0
      loff(i,2)=0
    1 continue
c
      do 2 i=me+1,np
      loff(i,1)=nqa
      loff(i,2)=nza
    2 continue
c
      call ga_igop(mrg_d01,loff,2*np,'+')
c
      ioffq=0
      if(me.gt.0) ioffq=loff(me,1)
      ioffz=loff(np,1)
      if(me.gt.0) ioffz=loff(np,1)+loff(me,2)
      nqat=loff(np,1)
      nzat=loff(np,2)
c
      return
      end
      subroutine qmmm_geometry(ndim,c,q,m,t,inum,ioffq,ioffz,
     + mwm,nwm,mwa,msa,nsa,iwz,isat,isdt,isq,isz,xw,xs)
c
      implicit none
c
#include "qmmm_common.fh"
#include "msgids.fh"
#include "global.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "constants.fh"
c
      integer cf_element
      real*8 cf_charge,cf_wcharge
      external cf_element,cf_charge,cf_wcharge
c
      integer ndim,mwm,nwm,mwa,msa,nsa,ioffq,ioffz
      real*8 c(3,ndim),q(ndim),m(ndim)
      character*16 t(ndim)
      integer inum(ndim)
      integer iwz(mwm),isat(msa),isdt(msa),isq(msa),isz(msa)
      real*8 xw(mwm,3,mwa),xs(msa,3)
c
      integer i,j,ioq,ioz
      real*8 atc,bqc
c
      ioq=ioffq
      ioz=ioffz
c
      do 1 i=1,ndim
      c(1,i)=0.0d0
      c(2,i)=0.0d0
      c(3,i)=0.0d0
      q(i)=0.0d0
      inum(i)=0
    1 continue
c
      atc=0.0d0
      bqc=0.0d0
      do 2 i=1,nsa
      if(iand(isdt(i),lquant).ne.0) then
      ioq=ioq+1
      c(1,ioq)=xs(i,1)*cnm2au
      c(2,ioq)=xs(i,2)*cnm2au
      c(3,ioq)=xs(i,3)*cnm2au
      inum(ioq)=cf_element(isat(i))
      q(ioq)=dble(inum(ioq))
      atc=atc+cf_charge(isq(i))
      endif
      if(isz(i).gt.0) then
      ioz=ioz+1
      c(1,ioz)=xs(i,1)*cnm2au
      c(2,ioz)=xs(i,2)*cnm2au
      c(3,ioz)=xs(i,3)*cnm2au
      if(iand(isdt(i),lqlink).ne.0) then
      q(ioz)=1.0d0
      inum(ioz)=1
      else
      q(ioz)=cf_charge(isq(i))
      bqc=bqc+q(ioz)
      inum(ioz)=0
      endif
      endif
    2 continue
      do 3 i=1,nwm
      if(iwz(i).gt.0) then
      do 4 j=1,mwa
      ioz=ioz+1
      c(1,ioz)=xw(i,1,j)*cnm2au
      c(2,ioz)=xw(i,2,j)*cnm2au
      c(3,ioz)=xw(i,3,j)*cnm2au
      q(ioz)=cf_wcharge(j)
      inum(ioz)=0
      bqc=bqc+q(ioz)
    4 continue
      endif
    3 continue
      bqc=bqc+dble(nint(atc))
c
      call ga_dgop(mrg_d02,c,3*ndim,'+')
      call ga_dgop(mrg_d03,q,ndim,'+')
      call ga_igop(mrg_d04,inum,ndim,'+')
      call ga_dgop(mrg_d05,bqc,1,'+')
      bqchrg=bqc
c
      if(.not.rtdb_put(irtdb,'charge',mt_dbl,1,bqchrg))
     + call errquit('Could not set charge',bqchrg)
c
      do 5 i=1,ndim
      call cf_num2tag(inum(i),t(i))
      if(inum(i).gt.0) then
      if(.not.geom_atn_to_default_mass(inum(i),m(i))) m(i)=0.0d0
      else
      m(i)=0.0d0
      endif
    5 continue
c
      return
      end
      subroutine qmmm_f(ndim,c,ioffq,ioffz,
     + mwm,nwm,mwa,msa,nsa,iwz,isat,isdt,isq,isz,fw,fs,energy)
c
      implicit none
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "qmmm_common.fh"
#include "constants.fh"
c
      integer ndim,mwm,nwm,mwa,msa,nsa,ioffq,ioffz
      real*8 c(3,ndim)
      integer iwz(mwm),isat(msa),isdt(msa),isq(msa),isz(msa)
      real*8 energy,fw(mwm,3,mwa),fs(msa,3)
c
      integer i,j,ioq,ioz
c
      ioq=ioffq
      ioz=ioffz
c
      if(.not.rtdb_get(irtdb,'task:gradient',mt_dbl,3*ndim,c))
     + call errquit('qmmm_f: rtdb_get gradient failed',0)
      if(.not.rtdb_get(irtdb,'task:energy',mt_dbl,1,energy))
     + call errquit('qmmm_f: rtdb_get energy failed',0)
c
c
      do 1 i=1,nsa
      if(iand(isdt(i),lquant).ne.0) then
      ioq=ioq+1
      fs(i,1)=fs(i,1)-c(1,ioq)*(cau2kj/cau2nm)
      fs(i,2)=fs(i,2)-c(2,ioq)*(cau2kj/cau2nm)
      fs(i,3)=fs(i,3)-c(3,ioq)*(cau2kj/cau2nm)
      endif
      if(isz(i).gt.0) then
      ioz=ioz+1
      fs(i,1)=fs(i,1)-c(1,ioz)*(cau2kj/cau2nm)
      fs(i,2)=fs(i,2)-c(2,ioz)*(cau2kj/cau2nm)
      fs(i,3)=fs(i,3)-c(3,ioz)*(cau2kj/cau2nm)
      endif
    1 continue
      do 2 i=1,nwm
      if(iwz(i).gt.0) then
      do 3 j=1,mwa
      ioz=ioz+1
      fw(i,1,j)=fw(i,1,j)-c(1,ioz)*(cau2kj/cau2nm)
      fw(i,2,j)=fw(i,2,j)-c(2,ioz)*(cau2kj/cau2nm)
      fw(i,3,j)=fw(i,3,j)-c(3,ioz)*(cau2kj/cau2nm)
    3 continue
      endif
    2 continue
c
      energy=energy*cau2kj
c
      return
      end
