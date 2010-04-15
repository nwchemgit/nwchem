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
c
      double precision energy
c
      logical task_fragment_energy_serial
      external task_fragment_energy_serial
c
      energy = 0.0d0
      pname = "task_fragment_energy"
c
      write(*,*) "initializing fragment"
      call smd_init("fragment",rtdb)
      call ga_sync()
      write(*,*) "finshed initializing fragment"
      return
c
      call util_message("task_smd_fragment_energy")
      result = task_fragment_energy_serial(rtdb)   

      if (.not. rtdb_put(rtdb,'fragment:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
      task_fragment_energy = .true.

      call util_message("done with task_smd_fragment_energy")
      return
      END

      function task_fragment_energy_serial(rtdb)
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
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision energy,ewald1,ewald2,evdw
      logical operiodic
      integer i
      integer group_id,group_nproc,ngroup
c
      call util_io_enable()
      call util_message("calling smd_monomer")
      call smd_monomer_energy(2,"dft")
      call smd_monomer_energy(1,"dft")
      call smd_monomer_energy(2,"dft")
c
      write(*,*) "Goodbye from task_smd"
      task_fragment_energy_serial = .true.
      return
      END
c
c      function task_smd_fragment_energy0(rtdb)
c      implicit none
c      integer rtdb
c
c#include "inp.fh"
c#include "mafdecls.fh"
c#include "rtdb.fh"
c#include "stdio.fh"
c#include "errquit.fh"
c#include "util.fh"
c#include "global.fh"
c      logical task_smd_fragment_energy0
c      character*72 namespace
c      character*32 pname
c      double precision etime1,etime2,etime3,etime4,etime5,etime6
c      double precision energy,ewald1,ewald2,evdw
c      logical operiodic
c      integer i
c      integer group_id,group_nproc,ngroup
c      double precision eps
cc
c      logical smd_tasklist_get1
c      external smd_tasklist_get1
c      logical smd_tasklist_get
c      external smd_tasklist_get
cc
c
c      pname = "task_smd_fragment_energy"
cc
c      eps = 1.0
cc
c      call smd_tasklist_init()
cc     
c      call ga_sync()
cc
c      call smd_group_init(rtdb)
cc
c      call smd_group_ngroup_get(ngroup)
c      call smd_group_id_get(group_id)
c      call smd_group_nproc_get(group_nproc)
c
c99    continue
c      if(.not.smd_tasklist_get1(i))
c     >  goto 100 
c
c      if(ga_nodeid().eq.0)
c     > write(*,*) "got task",i,group_id
c      call smd_monomer_energy(i)
cc
c      goto 99
c
c      call cpu_time(etime6)
c      write(luout,'(/,a,f20.3)')'Total CPU time : ',(etime6-etime1)
c      write(luout,'(/,a,f20.3)')'Initialization : ',(etime2-etime1)
c      write(luout,'(/,a,f20.3)')'Verlet         : ',(etime3-etime2)
c      write(luout,'(/,a,f20.3)')'Force          : ',(etime4-etime3)
c
c      write(*,*) "Goodbye from task_smd"
c100   continue
c      call ga_pgroup_sync(ga_pgroup_get_world())
c      call smd_fragment_charge_diff(eps)
c      if(eps.gt.0.0001) then
c        call smd_tasklist_set()
c        eps = 0.0d0
c        goto 99
c      end if
c
c      call ga_pgroup_sync(ga_pgroup_get_world())
c       call smd_group_end(rtdb)
c       energy = 0.0d0
c       if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
c     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
c      if(ga_nodeid().eq.0) then
c       write(*,*) "smd_charge_write_default()"
c       call smd_charge_write_default()
c      end if
c      call smd_end(rtdb)
c      task_smd_fragment_energy0 = .true.
c      return
c      END
c
c      function task_smd_fragment(rtdb)
c      implicit none
c      integer rtdb
c
c#include "inp.fh"
c#include "mafdecls.fh"
c#include "global.fh"
c#include "rtdb.fh"
c#include "stdio.fh"
c#include "errquit.fh"
c      logical task_smd_fragment
c      logical result,ofragment
cc
c      logical task_smd_fragment_energy
c      external task_smd_fragment_energy
c      logical task_smd_energy_sp
c      external task_smd_energy_sp
c      character*30 operation,pname
c      character*69 buffer
cc
c      double precision energy
cc
c      energy = 0.0d0
c      pname = "task_smd_fragment"
cc
c      write(*,*) "initializing fragment"
c      call smd_init("fragment",rtdb)
c      call ga_sync()
cc
c      if (.not. rtdb_cget(rtdb,'task:operation', 1,operation))
c     $      call errquit(pname,0, RTDB_ERR)
c
c      if(operation.eq."energy") then
c        call util_message("task_smd_fragment_energy")
c        result = task_smd_fragment_energy(rtdb)   
c      else
c        call errquit(pname//"unknown operation"//operation,
c     +               0, RTDB_ERR)
c      end if
c      if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
c     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
c      task_smd_fragment = .true.
c
c      return
c      END
c
c      function task_smd_fragment_energy(rtdb)
c      implicit none
c      integer rtdb
c
c#include "inp.fh"
c#include "mafdecls.fh"
c#include "rtdb.fh"
c#include "stdio.fh"
c#include "errquit.fh"
c#include "util.fh"
c#include "global.fh"
c      logical task_smd_fragment_energy
c      character*72 namespace
c      character*32 pname,rtdb_fname
c      double precision etime1,etime2,etime3,etime4,etime5,etime6
c      double precision energy,ewald1,ewald2,evdw
c      logical operiodic
c      integer i
c      integer group_id,group_nproc,ngroup
c      double precision eps
c      integer world_me
cc
c      logical smd_tasklist_get1
c      external smd_tasklist_get1
c      logical smd_tasklist_get
c      external smd_tasklist_get
cc
c
c      pname = "task_smd_fragment_energy"
cc
cc     closing world rtdb now
cc     its name is saved in smd data structure
cc     ---------------------------------------
c      call smd_rtdb_close('keep')
cc      if(.not.rtdb_close(rtdb, 'keep'))
cc     + call errquit('Failed to close world rtdb',0, RTDB_ERR)
c      call ga_sync() 
c      call smd_data_print(60+ga_nodeid())
c      call ga_sync() 
cc
cc     server-client loop
cc     ------------------
c      world_me = ga_nodeid()
c      
c      if(world_me.eq.0) then
c         call util_message("executing server code")
c         call sched_grp_server_code()
c      else
c         call util_message("executing client code")
c         call sched_grp_client_code()
c      endif
c
c      call ga_sync()
cc      call cpu_time(etime4)
c      call smd_rtdb_get_name(rtdb_fname)
c      if (.not. rtdb_open(rtdb_fname,'old',rtdb)) call errquit
c     *     (' rtdb_open old failed ',0,0)
cc
c
c
c      write(*,*) "ending smd"
c      call smd_data_print(30+ga_nodeid())
c      call ga_sync()
c      call smd_end(rtdb)
c      task_smd_fragment_energy = .true.
c      return
c      END
c
cc     ------------------------------
cc     simple subgroup implementation
cc     ------------------------------
c
c
