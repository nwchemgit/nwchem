c
c $Id: mm_interface.F,v 1.1 2004-02-10 01:07:32 marat Exp $
c
      subroutine mm_init(irtdb0)
      implicit none
c
#include "md_common.fh"
c
      integer irtdb0
c
c     set global rtdb for md module
c     ----------------------------
      irtdb=irtdb0
c
      call timer_init()
      call md_rdinp()
      call md_start()
c
      return
      end

      subroutine mm_end()
      implicit none
c
      call md_finish()
c
      return
      end
c
      subroutine mm_get_tot_nqm(nqm)
      implicit none
c
#include "mafdecls.fh"
#include "md_common.fh"
#include "msgids.fh"
c
      integer nqm
      integer tmp(1)
c
      call mm_get_nqm0(msa,nsaloc,
     &        int_mb(i_is+(lsdyn-1)*msa),
     &        tmp(1))
      call ga_igop(msg_qmmm_nat1,tmp,1,'+')
      nqm=tmp(1)
      end
c
      subroutine mm_get_nqm(nqm)
      implicit none
c
#include "mafdecls.fh"
#include "md_common.fh"
c
      integer nqm
c
      call mm_get_nqm0(msa,nsaloc,
     &        int_mb(i_is+(lsdyn-1)*msa),
     &        nqm)
      end
c
      subroutine mm_get_nqm0(msa,nsa,
     &                     isdt,nqm)
      implicit none
#include "util.fh"
#include "qmmm_common.fh"
c
      integer msa,nsa
      integer isdt(msa)
      integer nqm
c
      integer i
      integer nqatom,nlatom
c
      nqatom=0
      nlatom=0
      do i=1,nsa
        if(iand(isdt(i),lquant).ne.0) then
          nqatom=nqatom+1
        else if (and(isdt(i),lqlink).ne.0) then
          nlatom=nlatom+1
        endif
      end do
      nqm = nqatom + nlatom
c
c      if(util_print('qmmm',print_default)) then
c      write(*,8888) me,nlatom+nqatom
c      end if
c 8888 format(' Node ',1x,i1,'has ',i3,' quantum atoms')
c
      end
c
      subroutine mm_qm_mapping(nt,ai)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "md_common.fh"
      integer nt
      integer ai(nt)
c     local variables
      integer n
      integer l_ioff,i_ioff

      n=ga_nnodes()

      if(.not.ma_push_get(mt_int,n,'ioff',l_ioff,i_ioff))
     + call errquit('qmmm: Failed to allocate memory for ioff',np,
     &       MA_ERR)

      call mm_qm_mapping0(nt,n,msa,nsaloc,
     &                   int_mb(i_ioff),
     &                   int_mb(i_is+(lsdyn-1)*msa),
     &                   int_mb(i_is+(lsgan-1)*msa),
     &                   ai)

      if(.not.ma_pop_stack(l_ioff))
     + call errquit('qmmm: Failed to deallocate memory for ioff',np,
     &       MA_ERR)

      end
c
      subroutine mm_qm_mapping0(nt,np,msa,nsa,loff,isdt,isl,ai)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "qmmm_data.fh"
      integer nt,np
      integer msa,nsa
      integer loff(np)
      integer ai(nt)
      integer isl(msa)
      integer isdt(msa)

      integer i
      integer me
      integer nqa
      integer ioff
      integer ioq

      write(*,*) " mm_qm_mapping0"
      do i=1,np
        loff(i)=0
      end do
c
c     accumulate offsets
c     ------------------
      me = ga_nodeid()
      call mm_get_nqm(nqa)
      do i=me+1,np
        loff(i)=nqa
      end do    
      call ga_igop(msg_qmmm_offset,loff,np,'+')
c
c     set offset for a current node
c     -----------------------------
      if(me.gt.0) then
        ioff=loff(me)
      else
        ioff=0
      end if

      ioq=ioff
      do i=1,nsa
      if(iand(isdt(i),lquant).ne.0 .or.
     &   iand(isdt(i),lqlink).ne.0) then
         ioq=ioq+1
         ai(ioq)=isl(i)    
      endif
      end do
      do i=1,nt
        write(*,*) "i,ai",i,ai(i)
      end do
      end
c
c
      subroutine mm_get_qm_geom(nt,ai,inum,c,q,m,t)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "md_common.fh"
      integer nt
      integer ai(nt)
      integer inum(nt)
      double precision c(3,nt)
      double precision q(nt)
      double precision m(nt)
      character*16 t(nt)
c     local variables
      integer i

      call mm_get_qm_geom0(nt,msa,nsaloc,
     &                   int_mb(i_is+(lsgan-1)*msa),
     &                   int_mb(i_is+(lsdyn-1)*msa),
     &                   int_mb(i_is+(lsatt-1)*msa),
     &                   dbl_mb(i_xs),
     &                   ai,
     &                   inum,
     &                   c,
     &                   q,
     &                   m,
     &                   t)

      end
c
c
c
      subroutine mm_get_qm_geom0(nt,msa,nsa,
     &                   isl,
     &                   isdt,
     &                   isat,
     &                   xs,
     &                   ai,
     &                   inum,
     &                   c,
     &                   q,
     &                   m,
     &                   t)
      implicit none
#include "global.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "qmmm_data.fh"
#include "msgids.fh"
#include "inp.fh"
      integer nt,msa,nsa
      integer isl(msa)
      integer isdt(msa)
      integer isat(msa)
      double precision xs(msa,3)
      integer ai(nt)
      integer inum(nt)
      double precision c(3,nt)
      double precision q(nt)
      double precision m(nt)
      character*16 t(nt)
c
      integer cf_element
      external cf_element

c     local variables
      logical islink
      integer atn
      integer i,j
      character*16 tag

      write(*,*) "marat0"
      do i=1,nsa
         islink=iand(isdt(i),lquant).ne.0
         j=cf_element(isat(i))
         call cf_num2tag(j,tag)
         write(*,*) i,tag,xs(i,1),xs(i,2),xs(i,3),isl(i),islink
      end do

      write(*,*) "quantum atoms"
      call dfill(3*nt,0.0d0,c,1)
      do i=1,nsa
        islink=iand(isdt(i),lqlink).ne.0
        if(iand(isdt(i),lquant).ne.0 .or.
     &   islink) then
            do j=1,nt
              if(isl(i).eq.ai(j)) then
                 write(*,*) "found match",isl(i),ai(j)
                 write(*,*) "xs=",xs(i,1),xs(i,2),xs(i,3)
                 c(1,j)=xs(i,1)*cnm2au
                 c(2,j)=xs(i,2)*cnm2au
                 c(3,j)=xs(i,3)*cnm2au
                 atn=cf_element(isat(i))
                 if(islink) then
c                   if(atn.le.85) atn=85
c                   if(atn.le.53) atn=53
c                   if(atn.le.35) atn=35
c                   if(atn.le.17) atn=17
c                   if(atn.le.9)  atn=9
c                   if(atn.le.1)  atn=1
                    atn=1
                 end if
                 inum(j)=atn
              end if
            end do
        end if
      end do

      call ga_dgop(msg_qmmm_c,c,3*nt,'+')
      call ga_igop(msg_qmmm_inum,inum,nt,'+')

      do i=1,nt
         q(i) = dble(inum(i))
         call cf_num2tag(inum(i),t(i))
         if(.not.geom_tag_to_default_mass(t(i),m(i)))
     &   call errquit('qmmm: Failed to set masses',0, GEOM_ERR)
      end do

      end


cc
c      subroutine mm_get_qm_coord()
c      implicit none
cc
c#include "rtdb.fh"
c#include "util.fh"
c#include "mafdecls.fh"
c#include "md_common.fh"
cc
c      call qmmm_geom(irtdb,
c     & msa,nsaloc,int_mb(i_is+(lsatt-1)*msa),
c     & int_mb(i_is+(lsdyn-1)*msa),
c     & int_mb(i_is+(lsct1-1)*msa),int_mb(i_isz),
c     & dbl_mb(i_xs),dbl_mb(i_fs),
c     & linkatm)
cc
c      return
c      end
cc
cc
c      subroutine qmmm_geom(inrtdb
c     + msa,nsa,isat,isdt,isq,isz,xs,fs,energy,lnkatm)
c      implicit none
c#include "errquit.fh"
cc
c#include "rtdb.fh"
c#include "geom.fh"
c#include "util.fh"
c#include "global.fh"
c#include "mafdecls.fh"
c#include "qmmm_common.fh"
cc
c      external task_optimize
c      logical task_optimize
cc
c      integer cf_element
c      external cf_element
cc
c      integer inrtdb
c      integer mwm,nwm,mwa,nwa,msa,nsa,lnkatm
c      integer iwz(mwm),isat(msa),isdt(msa),isq(msa),isz(msa)
c      real*8 energy,xw(mwm,3,mwa),fw(mwm,3,mwa),xs(msa,3),fs(msa,3)
cc
c      integer i,np,ioffq,ioffz,igeom
c      integer i_ioff,l_ioff,i_c,l_c,i_q,l_q,i_t,l_t
c      integer i_i,l_i,i_j,l_j,i_m,l_m
c      integer nzatom,nqatom,nzatomt,nqatomt
c      logical xyz, xyzq
c      character*(nw_max_path_len) filename
cc
c      np=ga_nnodes()
c      me=ga_nodeid()
c      irtdb=inrtdb
cc
c      linkatm=lnkatm
cc
c      if(.not.ma_push_get(mt_int,2*np,'ioff',l_ioff,i_ioff))
c     + call errquit('qmmm: Failed to allocate memory for ioff',me,
c     &       MA_ERR)
cc
cc---- calculate number of quantum (nqatom) and Bq (nzatom) atoms.
c      nzatom=0
c      nqatom=0
c      do 1 i=1,nsa
c        if(iand(isdt(i),lquant).ne.0) then
c          nqatom=nqatom+1
c        endif
c    1 continue
cc
c      if(util_print('qmmm',print_default)) then
c      write(*,8888) me,nqatom
c      end if
c 8888 format(' Node ',1x,i1,'has ',i3,' quantum atoms')
c      stop
ccc
cc      call qmmm_offset(np,nqatom,nzatom,ioffq,ioffz,nqatomt,nzatomt,
cc     + int_mb(i_ioff))
ccc
cc      if(.not.ma_pop_stack(l_ioff))
cc     + call errquit('qmmm: Failed to deallocate memory for ioff',me,
cc     &       MA_ERR)
ccc
cc      if(.not.ma_push_get(mt_dbl,3*(nqatomt+nzatomt),'c',l_c,i_c))
cc     + call errquit('qmmm: Failed to allocate memory for c',me, MA_ERR)
cc      if(.not.ma_push_get(mt_dbl,nqatomt+nzatomt,'q',l_q,i_q))
cc     + call errquit('qmmm: Failed to allocate memory for q',me, MA_ERR)
cc      if(.not.ma_push_get(mt_dbl,nqatomt+nzatomt,'m',l_m,i_m))
cc     + call errquit('qmmm: Failed to allocate memory for m',me, MA_ERR)
cc      if(.not.ma_push_get(mt_byte,16*(nqatomt+nzatomt),'t',l_t,i_t))
cc     + call errquit('qmmm: Failed to allocate memory for t',me, MA_ERR)
cc      if(.not.ma_push_get(mt_int,nqatomt+nzatomt,'i',l_i,i_i))
cc     + call errquit('qmmm: Failed to allocate memory for i',me, MA_ERR)
cc      if(.not.ma_push_get(mt_int,2*(nqatomt+nzatomt),'j',l_j,i_j))
cc     + call errquit('qmmm: Failed to allocate memory for i',me, MA_ERR)
ccc
cc      call qmmm_geom_create(nqatomt+nzatomt,dbl_mb(i_c),dbl_mb(i_q),
cc     + dbl_mb(i_m),byte_mb(i_t),int_mb(i_i),int_mb(i_j),ioffq,ioffz,
cc     + mwm,nwm,mwa,msa,nsa,iwz,isat,isdt,isq,isz,xw,xs)
ccc
cc      if(.not.ma_pop_stack(l_j))
cc     + call errquit('qmmm: Failed to deallocate memory for j',me,
cc     &       MA_ERR)
cc      if(.not.ma_pop_stack(l_i))
cc     + call errquit('qmmm: Failed to deallocate memory for i',me,
cc     &       MA_ERR)
ccc
cc      if(.not.geom_create(igeom,'geometry'))
cc     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)
cc      if(.not.geom_cart_set(igeom,nqatomt+nzatomt,byte_mb(i_t),
cc     + dbl_mb(i_c),dbl_mb(i_q)))
cc     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)
cc      if(.not.geom_masses_set(igeom,nqatomt+nzatomt,dbl_mb(i_m)))
cc     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
cc      call geom_compute_values(igeom)
cc
cc      if (.not.rtdb_get(inrtdb,'qmmm:xyz',mt_log,1,xyz)) xyz=.false.
cc      if (.not.rtdb_get(inrtdb,'qmmm:xyzq',mt_log,1,xyzq)) xyzq=.false.
cc
ccc      if(xyz) then
cc        call util_file_name('xyz', .false., .false., filename)
cc        open(unit=43,status="new",form="formatted",file=filename)
cc        if(.not.geom_print_xyz(igeom,43))
cc     +    call errquit('qmmm: Failed to xyz file',0)
cc        close(43)
ccc      end if
cc
ccc      if(xyzq) then
cc        call util_file_name('xyzq', .false., .false., filename)
cc        open(unit=43,status="new",form="formatted",file=filename)
cc        if(.not.geom_print_xyzq(igeom,43))
cc     +    call errquit('qmmm: Failed to xyzq file',0)
cc        close(43)
ccc      end if
cc      stop
cc
ccc
cc      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
cc     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)
ccc
cc      if(.not.ma_pop_stack(l_t))
cc     + call errquit('qmmm: Failed to deallocate memory for t',me,
cc     &       MA_ERR)
cc      if(.not.ma_pop_stack(l_m))
cc     + call errquit('qmmm: Failed to deallocate memory for m',me,
cc     &       MA_ERR)
cc      if(.not.ma_pop_stack(l_q))
cc     + call errquit('qmmm: Failed to deallocate memory for q',me,
cc     &       MA_ERR)
ccc
cc      call ga_sync()
cc      if(.not.task_optimize(irtdb))
cc     + call errquit('qmmm: Task_gradient failed',0, GA_ERR)
cc      call ga_sync()
ccc
cc      call qmmm_f(nqatomt+nzatomt,dbl_mb(i_c),ioffq,ioffz,
cc     + mwm,nwm,mwa,msa,nsa,iwz,isat,isdt,isq,isz,fw,fs,energy)
ccc
cc      if(.not.ma_pop_stack(l_c))
cc     + call errquit('qmmm: Failed to deallocate memory for c',me,
cc     &       MA_ERR)
ccc
cc      if(.not.geom_destroy(igeom))
cc     + call errquit('qmmm: Failed to destroy geometry',me, GEOM_ERR)
ccc
cc      return
c      end
cc
c      subroutine qmmm_geom_create(ndim,c,q,m,t,inum,ndx,ioffq,ioffz,
c     + mwm,nwm,mwa,msa,nsa,iwz,isat,isdt,isq,isz,xw,xs)
cc
c      implicit none
c#include "errquit.fh"
cc
c#include "qmmm_common.fh"
c#include "msgids.fh"
c#include "global.fh"
c#include "geom.fh"
c#include "mafdecls.fh"
c#include "rtdb.fh"
c#include "constants.fh"
cc
c      integer cf_element
c      real*8 cf_charge,cf_wcharge
c      external cf_element,cf_charge,cf_wcharge
cc
c      integer ndim,mwm,nwm,mwa,msa,nsa,ioffq,ioffz
c      real*8 c(3,ndim),q(ndim),m(ndim)
c      character*16 t(ndim)
c      integer inum(ndim),ndx(ndim,2)
c      integer iwz(mwm),isat(msa),isdt(msa),isq(msa),isz(msa)
c      real*8 xw(mwm,3,mwa),xs(msa,3)
cc
c      integer i,j,ioq,ioz
c      real*8 atc,bqc,qqmatm
cc
c      ioq=ioffq
c      ioz=ioffz
cc
c      do 1 i=1,ndim
c      c(1,i)=0.0d0
c      c(2,i)=0.0d0
c      c(3,i)=0.0d0
c      q(i)=0.0d0
c      inum(i)=0
c      ndx(i,1)=0
c      ndx(i,2)=0
c    1 continue
cc
c      atc=0.0d0
c      bqc=0.0d0
c      do 2 i=1,nsa
c      if(iand(isdt(i),lquant).ne.0) then
c      ioq=ioq+1
c      c(1,ioq)=xs(i,1)*cnm2au
c      c(2,ioq)=xs(i,2)*cnm2au
c      c(3,ioq)=xs(i,3)*cnm2au
c      inum(ioq)=cf_element(isat(i))
c      ndx(ioq,1)=isdt(i)
c      ndx(ioq,2)=inum(ioq)
c      q(ioq)=dble(inum(ioq))
c      atc=atc+cf_charge(isq(i))
c      endif
c      if(isz(i).gt.0) then
c      ioz=ioz+1
c      c(1,ioz)=xs(i,1)*cnm2au
c      c(2,ioz)=xs(i,2)*cnm2au
c      c(3,ioz)=xs(i,3)*cnm2au
c      if(iand(isdt(i),lqlink).ne.0) then
c      if(linkatm.gt.0) then
c      q(ioz)=dble(linkatm)
c      inum(ioz)=linkatm
c      ndx(ioz,2)=linkatm
c      else
c      j=cf_element(isat(i))
c      ndx(ioz,2)=j
c      if(j.le.85) inum(ioz)=85
c      if(j.le.53) inum(ioz)=53
c      if(j.le.35) inum(ioz)=35
c      if(j.le.17) inum(ioz)=17
c      if(j.le.9) inum(ioz)=9
c      if(j.le.1) inum(ioz)=1
c      q(ioz)=dble(inum(ioz))
c      endif
c      else
c      q(ioz)=cf_charge(isq(i))
c      bqc=bqc+q(ioz)
c      inum(ioz)=0
c      endif
c      ndx(ioz,1)=isdt(i)
c      endif
c    2 continue
c      do 3 i=1,nwm
c      if(iwz(i).gt.0) then
c      do 4 j=1,mwa
c      ioz=ioz+1
c      c(1,ioz)=xw(i,1,j)*cnm2au
c      c(2,ioz)=xw(i,2,j)*cnm2au
c      c(3,ioz)=xw(i,3,j)*cnm2au
c      q(ioz)=cf_wcharge(j)
c      inum(ioz)=0
c      ndx(ioz,2)=0
c      bqc=bqc+q(ioz)
c    4 continue
c      endif
c    3 continue
cc
c      if(.not.rtdb_get(irtdb,'qmmm:qqmatm',mt_dbl,1,qqmatm)) then
c      bqc=bqc+dble(nint(atc))
c      else
c      if(me.eq.0) bqc=bqc+qqmatm
c      endif
cc
c      call ga_dgop(mrg_d02,c,3*ndim,'+')
c      call ga_dgop(mrg_d03,q,ndim,'+')
c      call ga_igop(mrg_d04,inum,ndim,'+')
c      call ga_igop(mrg_d05,ndx,2*ndim,'+')
c      call ga_dgop(mrg_d06,bqc,1,'+')
c      bqchrg=bqc
cc
c      if(.not.rtdb_put(irtdb,'charge',mt_dbl,1,bqchrg))
c     + call errquit('Could not set charge',bqchrg, RTDB_ERR)
cc
c      do 5 i=1,ndim
c      call cf_num2tag(inum(i),t(i))
c      if(iand(ndx(i,1),lqhigh).ne.0)
c     + t(i)=t(i)(1:index(t(i),' ')-1)//'_H'
c      if(iand(ndx(i,1),lqlink).ne.0)
c     + t(i)=t(i)(1:index(t(i),' ')-1)//'_L'
c      if(iand(ndx(i,1),lquant).ne.0.and.iand(ndx(i,1),lqsoft).ne.0)
c     + t(i)=t(i)(1:index(t(i),' ')-1)//'_S'
c      if(ndx(i,2).gt.0) then
c      if(.not.geom_atn_to_default_mass(ndx(i,2),m(i))) m(i)=0.0d0
c      else
c      m(i)=0.0d0
c      endif
c    5 continue
cc
c      return
c      end
c
