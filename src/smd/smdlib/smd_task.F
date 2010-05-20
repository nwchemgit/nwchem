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
      data itask /0/
      end
c
      subroutine smd_task_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_task
      character*32 tag,pname
      logical result

      pname = "smd_task_init_system"
c
      tag = "task"
      call util_message("in "//pname)
      call smd_system_get_component(sp_task,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_task_init(sp_task)

      return
      end
c      
      subroutine smd_task_init(namespace)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "smd_task_data.fh"
#include "smd_fids.fh"
#include "geom.fh"
c     
      character*(*) namespace
c
      character*32 pname
      character*80 tag
      character*80 tag1
      character*30 theory,operation
      integer nt,ns,i1,i2,i,j
      integer i_ifr,i_sop,i_sth,i_ifu
      logical result
      logical oparallel
c
      pname = "smd_task_init"
c
      call util_message("in "//pname)
c
c     create task namespace
c     ---------------------
      call smd_namespace_create(namespace)
c
c     find which monomers to do as given in the input file
c     ----------------------------------------------------
      tag1="smd:fragment:istart"
      call smd_rtdb_get_int(tag1,1,i1,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
      tag1="smd:fragment:iend"
      call smd_rtdb_get_int(tag1,1,i2,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
      nt = i2-i1+1
c
c     find theory and operation as given in task directive
c     ----------------------------------------------------
c      tag1="task:theory"
c     at this point theory has to be set explicitly
c     as it deafults to smd
c      tag1="smd:fragment:theory"
c      call smd_rtdb_get_string(tag1,1,theory,result)
c      if(.not.result) then
c        write(*,*) "need to set value for  "//tag1
c        call errquit("cannot find "//tag1,0,0)
c      end if

      tag1="task:operation"
      call smd_rtdb_get_string(tag1,1,operation,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
c
c     find if parallel mode is on
c     ---------------------------
      tag1="smd:fragment:parallel"
      call smd_rtdb_get_log(tag1,1,oparallel,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if

      tag1="task:operation"
      call smd_rtdb_get_string(tag1,1,operation,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
  
c
c     create task data structures
c     ---------------------------
      tag = "task:parallel"
      call smd_data_put_log(namespace,tag,1,oparallel)
      tag = "task:ntasks"
      call smd_data_put_int(namespace,tag,1,nt)
      tag = "task:func_id"
      call smd_data_create_get(namespace,tag,nt,MT_INT,i_ifu)
      tag = "task:frag_id"
      call smd_data_create_get(namespace,tag,nt,MT_INT,i_ifr)
c      tag = "task:theory"
c      call smd_data_create_get(namespace,tag,nt,MT_STR,i_sth)
      tag = "task:operation"
      call smd_data_create_get(namespace,tag,nt,MT_STR,i_sop)
c
c    fill up fragment id
c    -------------------
      j=0
      do i=i1,i2
cdebug
        int_mb(i_ifr+j) = i
        j=j+1
      enddo
c
      do i=1,nt
        int_mb(i_ifu+i-1) = i_smd_monomer
      enddo
c
c     fill in theory
c     --------------
c      call smd_set_char_to_ma2(theory,nt,byte_mb(i_sth))
c
c     fill in operation
c     -----------------
      call smd_set_char_to_ma2(operation,nt,byte_mb(i_sop))

      call util_message("out "//pname)

      return
      end
c
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
      character*30 pname
c
      pname = "smd_process_task"
      call util_message("in "//pname)
      if(ga_pgroup_get_default().ne.p_grp) 
     $  call errquit('Default group not properly set',p_grp,0)
c
      call cpu_time(etime1)
      call util_flush(6)
c      call smd_task_rtdb_clone(tskid,rtdb_grp)
      call util_flush(6)
      call cpu_time(etime2)
c      write(*,*) "time monomer rtdb",etime2-etime1

c      call smd_rtdb_set_handle(rtdb_grp)
c
      call ga_pgroup_sync(p_grp)
c      call smd_monomer(tskid)
      call ga_pgroup_sync(p_grp)
c
c      call smd_charge_ga_print()
      call ga_sync()
c      call smd_rtdb_close("delete")
c      if(.not.rtdb_close(rtdb_grp, 'delete'))
c     + call errquit('Failed to close group rtdb',0, GEOM_ERR)
c
      end subroutine
c
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
      call util_file_name_resolve(rtdb_name,.false.)
c     create name for group database in scratch
      write(buffer,'(I4.4,A3)') taskid,".db"
      call util_file_name(buffer,.true.,.false.,rtdb_name_group)
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
      return
      end
c
      subroutine smd_task_mode(oparallel)
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "smd_task_data.fh"
#include "smd_const_data.fh"

      logical oparallel
      character*(smd_string_size)  tag,sname
      logical result
      character*30 pname
c      
      pname = "smd_task_mode"
c      call util_message("in "//pname)
      tag = "task"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
c     get number of tasks
c     -------------------
      tag = "task:parallel"
      call smd_data_get_log(sname,tag,1,oparallel,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c      call util_message("out "//pname)
      return

      end subroutine
c
      subroutine smd_task_execute(itsk)
      implicit none
      integer itsk

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "smd_fids.fh"
#include "smd_const_data.fh"

      integer nt
      character*(smd_string_size)  tag,sname
      logical result
      character*30 pname
      integer i,nproc,ifunc
c      
      pname = "smd_task_execute"
      call util_message("in "//pname)
      tag = "task"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
c     get number of tasks
c     -------------------
      tag = "task:func_id"
      call smd_data_get_int1(sname,tag,itsk,itsk,ifunc)
      tag = "task:itsk"
      call smd_data_put_int(sname,tag,1,itsk)
      if(ifunc.eq.i_smd_monomer) then
       call smd_monomer()
      else
       call errquit(
     >       pname//'unknown function ',ifunc,0)
      end if
      call util_message("out "//pname)
      return

      end subroutine
c
      subroutine smd_task_iterate()
      implicit none
      integer itsk

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "smd_task_data.fh"
#include "smd_const_data.fh"
#include "msgids.fh"

      integer nt
      character*(smd_string_size)  tag,sname
      logical result
      character*30 pname
      integer i,j,nsb
      integer group_id
      integer ga_it
      logical odelay
      logical osb
c      
      pname = "smd_task_iterate"
c      call util_message("in "//pname)
      tag = "task"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
c     get number of tasks
c     -------------------
      tag = "task:ntasks"
      call smd_data_get_int(sname,tag,1,nt,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c      
      tag="smd:subgroups"
      call smd_rtdb_get_int(tag,1,nsb,osb)
c
      tag="smd:delay"
      call smd_rtdb_get_log(tag,1,odelay,result)
      if(.not.result) odelay=.false.
c
      if(.not.nga_create(mt_int, 1, 1, 
     >    'task counter',1, ga_it)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)


c     ---------------------------
      if(ga_nodeid().eq.0) then
         call nga_put(ga_it,1,1,1,1)
      end if
      call ga_sync()
 
      if(osb) then
      call smd_group_create_simple(nsb)
      call smd_group_id_get(group_id)
      if(ga_nodeid().eq.0) 
     +    write(luout,*) "Hello from group",group_id,i
      end if
c
10    continue    
      if(odelay) then
       if(group_id.ne.1) then
         write(luout,*) "sleeping in ",group_id
         call sleep(80)
       end if
      end if
      if(ga_nodeid().eq.0)
     >   j=NGA_Read_inc(ga_it, 1, 1)
      call ga_brdcst(msg_qmmm_misc,j,ma_sizeof(mt_int,1,mt_byte),0)
      if(j.gt.nt) goto 20
      if(ga_nodeid().eq.0) write(luout,*) "executing task",j
      call smd_task_execute(j)
      goto 10
c
20    continue
      if(osb) then
      call ga_pgroup_sync(ga_pgroup_get_world())
      call smd_group_end()
      end if
      if(.not.ga_destroy(ga_it)) call errquit(
     >       pname//'cannot destroy ga_it ',0,0)
c      call util_message("out "//pname)
      return

      end subroutine
c      
      subroutine smd_task_populate(task_list)
      implicit none
      integer task_list

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "smd_task_data.fh"
#include "smd_const_data.fh"

      integer nt
      character*(smd_string_size)  tag,sname
      logical result
      character*30 pname
      integer i,nproc
c      
      pname = "smd_task_populate"
      call util_message("in "//pname)
      tag = "task"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
c     get number of tasks
c     -------------------
      tag = "task:ntasks"
      call smd_data_get_int(sname,tag,1,nt,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

c      if(populate_done)
c     &  call errquit('trying to populate_tasks'// 
c     &    ' after done',0,0)

      nproc = 1
      if(populate_first_time) then

c       push in all available tasks
c       ---------------------------
        do i = 1,nt
c            call sched_grp_insert_task(task_list, i, nproc)
        enddo
        populate_first_time = .false.
      end if
      call util_message("out "//pname)
      return

      end subroutine
c
      subroutine smd_task_populate0(task_list)
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

      return
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
c            call sched_grp_insert_task(task_list, i, nproc)
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

