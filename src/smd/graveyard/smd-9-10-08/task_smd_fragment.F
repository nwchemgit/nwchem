      function task_smd_fragment(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      logical task_smd_fragment
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
      energy = 0.0d0
      pname = "task_smd_fragment"
c
      write(*,*) "initializing fragment"
      call smd_init("fragment",rtdb)
      call ga_sync()
c
      if (.not. rtdb_cget(rtdb,'task:operation', 1,operation))
     $      call errquit(pname,0, RTDB_ERR)

      if(operation.eq."energy") then
        call util_message("task_smd_fragment_energy")
        result = task_smd_fragment_energy(rtdb)   
      else
        call errquit(pname//"unknown operation"//operation,
     +               0, RTDB_ERR)
      end if
      if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
      task_smd_fragment = .true.

      return
      END

      function task_smd_fragment_energy(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
      logical task_smd_fragment_energy
      character*72 namespace
      character*32 pname,rtdb_fname
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision energy,ewald1,ewald2,evdw
      logical operiodic
      integer i
      integer group_id,group_nproc,ngroup
      double precision eps
      integer world_me
c
      logical smd_tasklist_get1
      external smd_tasklist_get1
      logical smd_tasklist_get
      external smd_tasklist_get
c

      pname = "task_smd_fragment_energy"
c
c     closing world rtdb now
c     its name is saved in smd data structure
c     ---------------------------------------
      call smd_rtdb_close('keep')
c      if(.not.rtdb_close(rtdb, 'keep'))
c     + call errquit('Failed to close world rtdb',0, RTDB_ERR)
      call ga_sync() 
      call smd_data_print(60+ga_nodeid())
      call ga_sync() 
c
c     server-client loop
c     ------------------
      world_me = ga_nodeid()
      
      if(world_me.eq.0) then
         call util_message("executing server code")
         call sched_grp_server_code()
      else
         call util_message("executing client code")
         call sched_grp_client_code()
      endif

      call ga_sync()
c      call cpu_time(etime4)
      call smd_rtdb_get_name(rtdb_fname)
      if (.not. rtdb_open(rtdb_fname,'old',rtdb)) call errquit
     *     (' rtdb_open old failed ',0,0)
c


      write(*,*) "ending smd"
      call smd_data_print(30+ga_nodeid())
      call ga_sync()
      call smd_end(rtdb)
      task_smd_fragment_energy = .true.
      return
      END

c     ------------------------------
c     simple subgroup implementation
c     ------------------------------
      function task_smd_fragment_energy0(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
      logical task_smd_fragment_energy0
      character*72 namespace
      character*32 pname
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision energy,ewald1,ewald2,evdw
      logical operiodic
      integer i
      integer group_id,group_nproc,ngroup
      double precision eps
c
      logical smd_tasklist_get1
      external smd_tasklist_get1
      logical smd_tasklist_get
      external smd_tasklist_get
c

      pname = "task_smd_fragment_energy"
c
      eps = 1.0
c
      call smd_tasklist_init()
c     
      call ga_sync()
c
      call smd_group_init(rtdb)
c
      call smd_group_ngroup_get(ngroup)
      call smd_group_id_get(group_id)
      call smd_group_nproc_get(group_nproc)

99    continue
      if(.not.smd_tasklist_get1(i))
     >  goto 100 

      if(ga_nodeid().eq.0)
     > write(*,*) "got task",i,group_id
      call smd_monomer_energy(i)
c
      goto 99

      call cpu_time(etime6)
      write(luout,'(/,a,f20.3)')'Total CPU time : ',(etime6-etime1)
      write(luout,'(/,a,f20.3)')'Initialization : ',(etime2-etime1)
      write(luout,'(/,a,f20.3)')'Verlet         : ',(etime3-etime2)
      write(luout,'(/,a,f20.3)')'Force          : ',(etime4-etime3)

      write(*,*) "Goodbye from task_smd"
100   continue
      call ga_pgroup_sync(ga_pgroup_get_world())
      call smd_fragment_charge_diff(eps)
      if(eps.gt.0.0001) then
        call smd_tasklist_set()
        eps = 0.0d0
        goto 99
      end if

      call ga_pgroup_sync(ga_pgroup_get_world())
       call smd_group_end(rtdb)
       energy = 0.0d0
       if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
      if(ga_nodeid().eq.0) then
       write(*,*) "smd_charge_write_default()"
       call smd_charge_write_default()
      end if
      call smd_end(rtdb)
      task_smd_fragment_energy0 = .true.
      return
      END


c $Id$
