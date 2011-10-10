      block data process_data
      implicit none
#include "smd_task_data.fh"
      data populate_first_time /.true./
      data populate_done /.false./
      data itask /0/
      end block data
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

c      call smd_task_init(sp_task)

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
c      tag1="smd:fragment:parallel"
c      call smd_rtdb_get_log(tag1,1,oparallel,result)
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
c     create task data structures
c     ---------------------------
c      tag = "task:parallel"
c      call smd_data_put_log(namespace,tag,1,oparallel)
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
       call smd_task_init1()
       call smd_monomer()
       call smd_task_end1()
      else
       call errquit(
     >       pname//'unknown function ',ifunc,0)
      end if
      call util_message("out "//pname)
      return

      end subroutine
c
      subroutine smd_task_get_counter(im)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "smd_const_data.fh"
      integer im
c
      character*32 pname
      logical oscratch
      character*(smd_string_size)  tag,sname
      logical result 
      integer itsk

      pname = "smd_monomer_init"
c
      tag = "task"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
c     get actual task reference
c     -------------------------
      tag = "task:itsk"
      call smd_data_get_int(sname,tag,1,itsk,result)
      if(.not.result)
     *     call errquit('cannot find'//tag,0,0)
c      write(luout,*) "found itsk",itsk
c      
      tag = "task:frag_id"
      call smd_data_get_int1(sname,tag,itsk,itsk,im)
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
      integer iskip
c      
      call util_io_enable()
      pname = "smd_task_iterate"
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
c      
      tag="smd:subgroups"
      call smd_rtdb_get_int(tag,1,nsb,osb)
c
      tag="smd:delay"
      call smd_rtdb_get_log(tag,1,odelay,result)
      if(.not.result) odelay=.false.
c
      tag="smd:iskip"
      call smd_rtdb_get_int(tag,1,iskip,result)
      if(.not.result) iskip=0
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
      if(group_id.eq.iskip) goto 20
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
      if(ga_nodeid().eq.0) write(luout,*) "executing task",j,group_id
      call smd_task_execute(j)
      goto 10
c
20    continue
      call util_message("1 "//pname)
      if(osb) then
      call ga_pgroup_sync(ga_pgroup_get_world())
      call smd_group_end()
      end if
      if(.not.ga_destroy(ga_it)) call errquit(
     >       pname//'cannot destroy ga_it ',0,0)
      call util_message("out "//pname)
      return

      end subroutine
c      
      subroutine smd_task_init1()
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
c
      character*32 pname
      logical oscratch
      character*(smd_string_size)  tag
      character*(smd_string_size)  sp_geom
      character*(smd_string_size)  sp_task
      character*(smd_string_size)  operation
      character*(smd_string_size)  theory
      logical result 
      character*64 msg
      integer im
      logical oinclude
      integer rtdb

      pname = "smd_task_init1"
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "task"
      call smd_system_get_component(sp_task,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c      
c     current task index
c     ------------------
      tag = "task:itsk"
      call smd_data_get_int(sp_task,tag,1,itsk,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no task:itsk '//tag,0,0)
c
c     grab global charges and coords
c     ------------------------------
c      call smd_global_update("charges","get")
c      call smd_global_update("coords","get")
c     set central residue
      tag = "task:frag_id"
      call smd_data_get_int1(sp_task,tag,itsk,itsk,im)
      if(ga_nodeid().eq.0) then
        write(*,*) "preparing for fragment",im
      end if
      tag = "icenter"
      call smd_data_put_int(sp_geom,tag,1,im)

      tag = "oinclude"
      call smd_data_get_log2(sp_task,tag,itsk,itsk,oinclude,result)
      if(.not.result) oinclude=.true.
      call smd_data_put_log(sp_geom,tag,1,oinclude)
c
c     create a copy of rtdb
c     ---------------------
      oscratch = .true.
      call smd_rtdb_push_open_clone(itsk,oscratch,rtdb)
c
c     set operation theory
c     --------------------
      tag="task:operation"
      call smd_data_get_string2(sp_task,tag,itsk,itsk,
     >                          operation,result)
      if(result) then
        if (.not.(rtdb_cput(rtdb,'task:operation',1,operation)))
     *       call errquit('rtdb_cput task:operation failed',0,0)
      end if

      tag="task:theory"
      call smd_data_get_string2(sp_task,tag,itsk,itsk,
     >                          theory,result)
      if(result) then
        if (.not.(rtdb_cput(rtdb,'task:theory',1,theory)))
     *       call errquit('rtdb_cput task:theory failed',0,0)
      end if

c     create geometry object
c     ---------------------
      call smd_geom_prepare()
      call smd_bq_init()
      call ga_sync()
      call smd_group_set_io_custom(im,"frag")
      return
      end subroutine
C>    Postprocessing and cleanup after
C>    task has been performed
C>    @see smd_task_init1() for preprocessing
      subroutine smd_task_end1()
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
#include "bq.fh"
c
      character*32 pname
      logical oscratch
      character*(smd_string_size)  tag
      character*(smd_string_size)  sp_geom
      character*(smd_string_size)  sp_task
      character*(smd_string_size)  operation
      character*(smd_string_size)  theory
      logical result 
      character*64 msg
      integer im
      logical oinclude
      integer rtdb
      integer bq_handle

      pname = "smd_task_end1"
c
c      write(*,*) pname//"10"
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "task"
      call smd_system_get_component(sp_task,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c      
c     current task index
c     ------------------
      tag = "task:itsk"
      call smd_data_get_int(sp_task,tag,1,itsk,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no task:itsk '//tag,0,0)
c     central residue
      tag = "task:frag_id"
      call smd_data_get_int1(sp_task,tag,itsk,itsk,im)
c
c     what was the operation
c     ----------------------
      call smd_rtdb_get_handle(rtdb)
      if (.not.(rtdb_cget(rtdb,'task:operation',1,operation)))
     *       call errquit('rtdb_cget task:operation failed',0,0)
c      
      if(inp_compare(.false., 'dipole', operation)) then
        call smd_geom_get_dipole()
        call smd_global_update_res("dipoles","put",im)
      else if(inp_compare(.false., 'gradient', operation)) then
        call smd_geom_get_gradients()
        call smd_global_update_res("gradients","put",im)
      else if(inp_compare(.false., 'optimize', operation)) then
        call smd_geom_get_coords()
        call smd_global_update_res("coords","put",im)
      else if(inp_compare(.false., 'espcharge', operation)) then
        if(ga_nodeid().eq.0) then
          write(luout,*) "updating charges for fragment",im
        end if
        call smd_geom_get_charges()
      end if

      call ga_sync()
      if(bq_get_active(bq_handle)) then
        if(.not.bq_destroy(bq_handle))
     +   call errquit(pname//'Failed bq_destroy',0,CALC_ERR)
      end if

      call smd_geom_end()
      call smd_rtdb_pop()

c      write(*,*) pname//"100"
 
      call smd_group_unset_io()
      return
      end subroutine
c

c $Id$
