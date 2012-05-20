c     Wrapper method to invoke calculation on an SMD monomer calculation.
c     @author Sriram Krishnamoorthy
c
c     This method is invoked from the dynamic process group client
C     -server implementation.  


c     Wrapper method invoked from the C client-server scheduler
c     @param tskid IN Id of task to execute (here monomer id)
c     @param p_grp IN Handle to GA group on which to execute
c

      block data process_data
      implicit none
#include "smd_task_data.fh"
      data populate_first_time /.true./
      data populate_done /.false./
      end

      subroutine smd_process_task(tskid, p_grp)
c
c     commences task(taskid) for a group (p_grp)
      implicit none
      integer tskid, p_grp 

#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"

      integer rtdb_grp 
      double precision etime1,etime2
c
      if(ga_pgroup_get_default().ne.p_grp) 
     $  call errquit('Default group not properly set',p_grp,0)
c
      call cpu_time(etime1)
      write(*,*) "started rtdb cloning"
      call util_flush(6)
      call smd_task_rtdb_clone(tskid,rtdb_grp)
      write(*,*) "finished rtdb cloning"
      call util_flush(6)
      call cpu_time(etime2)
c      write(*,*) "time monomer rtdb",etime2-etime1

c      call smd_rtdb_set_handle(rtdb_grp)
c
      call ga_pgroup_sync(p_grp)
      call smd_monomer_energy(tskid)
      call ga_pgroup_sync(p_grp)
c
c      call smd_charge_ga_print()
      call ga_sync()
      call smd_rtdb_close("delete")
c      if(.not.rtdb_close(rtdb_grp, 'delete'))
c     + call errquit('Failed to close group rtdb',0, GEOM_ERR)
c
      end subroutine

      subroutine smd_task_rtdb_clone(taskid,rtdb)
c     rtdb cloning
c     1. make a numbered name (taskid)
c     2. copy file
c     3. open rtdb
c
c     taskid[IN] task number
c     rtdb[OUT]  rtdb handle
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
#include "smd_const_data.fh"
#include "rtdb.fh"
      integer taskid,rtdb
c     
      character*(smd_string_size)  buffer
      character*256  rtdb_name
      character*256  rtdb_name_group
      integer i1,i2
      logical result
      logical master
      character*30 pname
c
      pname = "smd_task_rtdb_clone"
c
      call ga_sync()
      master = ga_nodeid().eq.0
c
c     geti absolute name of the original database
      call smd_rtdb_get_name(rtdb_name)
      write(*,*) "rtdb_name",rtdb_name
      call util_file_name_resolve(rtdb_name,.false.)
      write(*,*) "resolved rtdb_name",rtdb_name
c     create name for group database in scratch
      write(buffer,'(I4.4,A3)') taskid,".db"
      call util_file_name(buffer,.true.,.false.,rtdb_name_group)
      write(*,*) "resolved group_rtdb_name",rtdb_name_group
      call util_flush(6)
c
      if(master) then
        inquire(file=rtdb_name,exist=result)
        if(.not.result) then
          call errquit("cannot find rtdb"//
     >                 rtdb_name,
     >                 0,0)
        end if
        call util_file_copy(rtdb_name,rtdb_name_group)
      end if

      call smd_rtdb_open(rtdb_name_group,"old")
c      if(.not.rtdb_open(rtdb_name_group,"old",rtdb))
c     + call errquit('Failed to open group rtdb',0, GEOM_ERR)
      write(*,*) "got new rtdb handle",rtdb

      return
      end


c     task_list IN Task list to insert tasks into
      subroutine smd_task_populate(task_list)
      integer task_list

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "smd_task_data.fh"

      integer nr, na, i, ntsks, ga_charges
      logical changed_flag
      integer lo(1), hi(1), ld(1)
      double precision diff
      integer istart,iend,nrpoc

c      istart=1
c      iend=4
      call smd_rtdb_get_istart(istart)
      call smd_rtdb_get_iend(iend)
      call smd_rtdb_get_nproc(nproc)
      if(populate_done)
     &  call errquit('trying to populate_tasks'// 
     &    ' after done',0,0)

      call smd_system_ntot(na)
      call smd_system_nres(nr)


      ntsks = 0
      if(populate_first_time) then

c         call smd_charge_get_ptr(i_q)
c
c         if(.not.ma_alloc_get(MT_DBL,
c     +        na,
c     +        "tmp q",
c     +        tmp_h_q,
c     +        tmp_i_q))
c     +        call errquit('memory allocation',
c     +        na,
c     +        MA_ERR)
c
c         old_i_q = i_q
c         new_i_q = tmp_i_q

c         do i = 1, nr
        do i = istart,iend
            call sched_grp_insert_task(task_list, i, nproc)
            ntsks = ntsks + 1
         enddo
         
         populate_first_time = .false.
      else
c
c         call smd_charge_ga_handle(ga_charges)
c
c         lo(1) = 1
c         hi(1) = na
c         ld(1) = 1
c         call nga_get(ga_charges,lo, hi, dbl_mb(new_i_q), ld)
c
c         changed_flag = .false.
c         do i = 1, na
c            diff = dabs(dbl_mb(new_i_q+i-1)-dbl_mb(old_i_q+i-1))
c            if(diff.gt.0.00001) then
c              write(*,*) "charge difference",
c     >                   diff
c              changed_flag = .true.
c              goto 100
c            end if
c         end do
c 100     continue
c
c         if(changed_flag) then
cc            do i = 1, nr
c           do i = istart,iend
c               call sched_grp_insert_task(task_list, i, 1)
c               ntsks = ntsks + 1
c            enddo
c         endif
c
c         !Flip new and old pointers
c         i = new_i_q
c         new_i_q = old_i_q
c         old_i_q = i
      end if

c      if(ntsks.eq.0) then
c         if(.not.ma_free_heap(tmp_h_q))
c     $        call errquit('memory deallocation',0,tmp_h_q)
c         populate_done=.true.
c      end if

      end !populate_tasks

c $Id$
