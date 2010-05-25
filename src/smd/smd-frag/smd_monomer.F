      function smd_monomer()
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
      logical smd_monomer
      integer im
      character*(30) theory
c
      integer rtdb
      character*72 operation
      character*32 pname
      logical result
      logical task_energy_doit
      external task_energy_doit
      logical task_optimize
      external task_optimize
      logical task_gradient
      external task_gradient
      logical task_energy
      external task_energy
      logical esp
      external esp
      integer util_sgroup_mygroup
      external util_sgroup_mygroup
      character*255 key
      logical noesp

      pname = "smd_monomer"

      call smd_monomer_init(im,theory,operation,rtdb)
      if(inp_compare(.false., 'energy', operation)) then
         result = task_energy(rtdb)
      else if(inp_compare(.false., 'gradient', operation)) then
        result = task_gradient(rtdb)
        call smd_geom_get_gradients()
        call smd_geom_get_coords()
        call smd_global_update_res("gradients","put",im)
        call ga_sync()
c        call smd_global_print("gradients")
      else if(inp_compare(.false., 'optimize', operation)) then
        result = task_optimize(rtdb)
        call smd_geom_get_coords()
        call smd_global_update_res("coords","put",im)
        call ga_sync()
c        call smd_global_print("coords")
      else if(inp_compare(.false., 'espcharge', operation)) then
         result = task_energy(rtdb)
         if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
     +        call errquit('esp: rtdb_put newgrd failed',0,0)
         if(.not.esp(rtdb))
     &    call errquit(pname//'failed esp call', 100,
     &                 RTDB_ERR)
          call smd_geom_get_charges()
          if(ga_nodeid().eq.0)
     &         call smd_fragment_charge_set(im)
      else
        call errquit(pname//"unknown operation "//operation,0,0)
      end if
      call smd_coords_save(im,.true.)
cc
c      if(.not.rtdb_get(rtdb,'smd:noesp',mt_log,1,noesp))
c     +     noesp = .false.
cc
c      if(noesp) goto 100
cc
c      if(theory.ne."dft") goto 100
c      if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
c     +     call errquit('esp: rtdb_put newgrd failed',0,0)
c
c
c      if(.not.esp(rtdb))
c     &    call errquit(pname//'failed esp call', 100,
c     &                 RTDB_ERR)
cc
c      call util_message("finished esp")
cc
c      call smd_geom_get_charges()
c      if(ga_nodeid().eq.0)
c     &     call smd_fragment_charge_set(im)
c 
c      call util_message("finished smd_geom_get_charges")
c
100   continue

      write(key,'(a,'':converged'')') theory(1:inp_strlen(theory))
      if (.not. rtdb_put(rtdb, key, mt_log, 1, .false.))
     &    call errquit(pname//'  rtdb_put', 0,
     &                 RTDB_ERR)
c      call smd_fragment_charge_print()

      smd_monomer = .true.

      call smd_monomer_end()
      call util_message("finished "//pname)
      return
      end function

      subroutine smd_monomer_init(im,theory,operation,rtdb)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "smd_const_data.fh"
      integer im
      character*(*) theory
      character*(*) operation
      integer rtdb
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
      write(luout,*) "found itsk",itsk
c      
      tag="task:operation"
      call smd_data_get_string1(sname,tag,itsk,itsk,
     >                          operation)
c
      tag = "task:frag_id"
      call smd_data_get_int1(sname,tag,itsk,itsk,im)
c
      tag = "task:itsk"
      call smd_data_get_int(sname,tag,1,itsk,result)
      if(.not.result)
     *     call errquit('cannot find'//tag,0,0)
c
      oscratch = .true.
      call smd_rtdb_push_open_clone(itsk,oscratch,rtdb)
c     get the theory 
c     --------------
      if (.not.(rtdb_cget(rtdb,'task:theory',1,theory)))
     *     call errquit('rtdb_cget task:theory failed',0,0)
c
c     pull the global charges/coords onto local set
      call smd_fragment_charge_get_all()
      call smd_global_update("coords","get")
c
c      create geometry object
c      ---------------------
      call smd_geom_init_system(im)
      call smd_bq_init()
      call ga_sync()
c
      call smd_message_post("out "//pname)

      return
      end subroutine

      subroutine smd_monomer_end()
      implicit none
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "bq.fh"
#include "smd_const_data.fh"

      character*30 pname
      integer bq_handle
    
      pname = "smd_monomer_end"
      if(bq_get_active(bq_handle)) then
        if(.not.bq_destroy(bq_handle))
     +   call errquit(pname//'Failed bq_destroy',0,CALC_ERR)
      end if

      call smd_rtdb_pop()
 
      return
      end subroutine

