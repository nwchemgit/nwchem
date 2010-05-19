      function smd_monomer(itsk)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
      logical smd_monomer
      integer im,itsk
      character*(30) theory
c
      integer rtdb
      character*72 operation
      character*72 namespace
      character*32 pname
      logical result
      double precision qm_energy
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
      logical mode
      character*72 tag
      character*256 movecs
      character*256 rtdb_name
      double precision e1, e2,e3,e4,e5,e6,e7,e8,e9,e10,e11
      character*255 prefix, prefix0,key
      integer i1
      integer group_id
      logical noesp

      pname = "smd_monomer"

      call smd_monomer_init(itsk,im,theory,operation,rtdb)
c      write(*,*) "operation is", operation
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
       else
        call errquit(pname//"unknown operation "//operation,0,0)
      end if
c
c      if(ga_nodeid().eq.0)
c     >  write(*,*) "finished energy",theory,qm_energy
c      
      if(.not.rtdb_get(rtdb,'smd:noesp',mt_log,1,noesp))
     +     noesp = .false.
c
      if(noesp) goto 100
c
      if(theory.ne."dft") goto 100
      if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
     +     call errquit('esp: rtdb_put newgrd failed',0,0)


      if(.not.esp(rtdb))
     &    call errquit(pname//'failed esp call', 100,
     &                 RTDB_ERR)
c
      call util_message("finished esp")
c
      call smd_geom_get_charges()
      if(ga_nodeid().eq.0)
     &     call smd_fragment_charge_set(im)
 
      call util_message("finished smd_geom_get_charges")

100   continue
      call util_message("finished smd_geom_monomer_end")

      write(key,'(a,'':converged'')') theory(1:inp_strlen(theory))
      if (.not. rtdb_put(rtdb, key, mt_log, 1, .false.))
     &    call errquit(pname//'  rtdb_put', 0,
     &                 RTDB_ERR)
      call smd_fragment_charge_print()

      smd_monomer = .true.

      call smd_geom_monomer_end()
      call util_message("finished "//pname)
      return
      end function

      subroutine smd_monomer_init(itsk,im,theory,operation,rtdb)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "smd_const_data.fh"
      integer itsk
      integer im
      character*(*) theory
      character*(*) operation
      integer rtdb
c
      character*72 namespace
      character*32 pname
      integer group_id
      logical oscratch
      character*(smd_string_size)  tag,sname,buffer
      logical result 

      pname = "smd_monomer_init"
c
      tag = "task"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c      
      tag="task:operation"
      call smd_data_get_string1(sname,tag,itsk,itsk,
     >                          operation)
c      write(*,*) "found operation",operation
c
c      tag="task:theory"
c      call smd_data_get_string1(sname,tag,itsk,itsk,
c     >                          theory)
c      write(*,*) "found theory",theory
c
      tag = "task:frag_id"
      call smd_data_get_int1(sname,tag,itsk,itsk,im)
c
      oscratch = .true.
      call smd_rtdb_push_open_clone(itsk,oscratch,rtdb)
      write(luout,*) "now all my output hsould go here" 
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
c      call smd_geom_monomer_create(im)
      call ga_sync()
c
      call smd_message_post("out "//pname)

      return
      end subroutine


