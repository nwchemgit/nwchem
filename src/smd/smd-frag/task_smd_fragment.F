      function task_fragment_energy(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      logical task_fragment_energy
      logical result,ofragment
c
      logical task_smd_fragment_energy
      external task_smd_fragment_energy
      logical task_smd_energy_sp
      external task_smd_energy_sp
      character*30 operation,pname
      character*69 buffer
      integer world_me
c
      double precision energy
      logical oparallel
c
      logical  task_fragment_energy_parallel
      external task_fragment_energy_parallel
      logical task_fragment_energy_serial
      external task_fragment_energy_serial
c
      pname = "task_fragment_energy"
c
      call smd_init("fragment",rtdb)
      call ga_sync()
c
c     close main rtdb as we may be going into groups
c     ----------------------------------------------
      call smd_rtdb_close("keep")
c
      energy = 0.0d0
      if(ga_nodeid().eq.0) 
     +    call util_io_enable()
      call smd_task_mode(oparallel)
      if (oparallel) then
        call util_message("executing in parallel mode")
        result = task_fragment_energy_parallel()
      else
        call util_message("executing in serial mode")
        result = .true.
        call smd_task_iterate()
c        result = task_fragment_energy_serial()   
      end if
      call smd_rtdb_open(rtdb)
      energy =0.0d0
      if (.not. rtdb_put(rtdb,'fragment:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
      task_fragment_energy = .true.

      call util_message("done with task_smd_fragment_energy")
      if(ga_nodeid().eq.0) then
       call smd_charge_write_default()
      end if
      call smd_global_end()
      call smd_end()
      return
      END

      function task_fragment_energy_parallel()
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      logical task_fragment_energy_parallel
      logical result,ofragment
c
      character*30 operation,pname
      character*69 buffer
      integer world_me
c
c
c
      pname = "task_fragment_energy_parallel"
c     testing
c     server-client loop
c     ------------------
      task_fragment_energy_parallel =.true.
      world_me = ga_nodeid()
      
      if(world_me.eq.0) then
         call util_message("executing server code")
c         call sched_grp_server_code()
      else
         call util_message("executing client code")
c         call sched_grp_client_code()
      endif

      call ga_sync()
      return
c     end testing
      call util_message("out "//pname)
      return
      END

      function task_fragment_energy_serial()
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
      logical task_fragment_energy_serial
      character*72 namespace
      character*32 pname
c
      call util_io_enable()
      call util_message("calling smd_monomer")
      call smd_monomer(2)
      call smd_monomer(1)
      call smd_monomer(2)
c
      write(*,*) "Goodbye from task_smd"
      task_fragment_energy_serial = .true.
      return
      END

      subroutine task_fragment_populate()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "smd_fids.fh"
c     
      character*32 pname
      character*80 tag
      character*80 namespace
      character*80 tag1
      character*30 operation
      integer nt,ns,i1,i2,i,j
      integer i_ifr,i_sop,i_sth,i_ifu
      logical result
      logical oparallel
c
      pname = "task_fragment_populate"
c
      call util_message("in "//pname)
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
      tag1="task:operation"
      call smd_rtdb_get_string(tag1,1,operation,result)
      if(.not.result) then
        write(*,*) "need to set value for  "//tag1
        call errquit("cannot find "//tag1,0,0)
      end if
c
c     create task namespace
c     ---------------------
      namespace = "fragment task"
      call smd_namespace_create(namespace)
c
c     create task data structures
c     ---------------------------
      tag = "task:ntasks"
      call smd_data_put_int(namespace,tag,1,nt)
      tag = "task:func_id"
      call smd_data_create_get(namespace,tag,nt,MT_INT,i_ifu)
      tag = "task:frag_id"
      call smd_data_create_get(namespace,tag,nt,MT_INT,i_ifr)
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
c     fill in operation
c     -----------------
      call smd_set_char_to_ma2(operation,nt,byte_mb(i_sop))

      call util_message("out "//pname)

      return
      end
