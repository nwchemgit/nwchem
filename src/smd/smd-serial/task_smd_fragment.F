      function task_smd_fragment(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      logical task_smd_fragment
      logical ogroups,result,ofragment
c
      logical task_smd_fragment_energy
      external task_smd_fragment_energy
      logical task_smd_energy_sp
      external task_smd_energy_sp
      character*30 operation,pname
c
      double precision energy
c
      energy = 0.0d0
      pname = "task_smd_fragment"
c
      call smd_groups_inquire(rtdb,ogroups)
c
      call smd_init("fragment",rtdb)
c
      if (.not. rtdb_cget(rtdb,'task:operation', 1,operation))
     $      call errquit(pname,0, RTDB_ERR)

      write(*,*) "found operation",operation
      if(operation.eq."energy") then
        result = task_smd_fragment_energy(rtdb)   
      else
        call errquit(pname//"unknown operation"//operation,
     +               0, RTDB_ERR)
      end if
       write(*,*) "Goodbye from task_smd_fragment"
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
      character*32 pname
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision energy,ewald1,ewald2,evdw
      logical operiodic
      integer i
c
      integer util_sgroup_mygroup
      external util_sgroup_mygroup

      logical smd_tasklist_get
      external smd_tasklist_get
c
      integer util_sgroup_numgroups
      external util_sgroup_numgroups

      pname = "task_smd_fragment_energy"
c
      call smd_tasklist_init()
c     
      call ga_sync()
c
      call smd_init_groups(rtdb)
c
      write(*,*) "number of groups",util_sgroup_numgroups(),
     >           util_sgroup_mygroup()

99    continue
      if(.not.smd_tasklist_get(i))
     >  goto 100 
      if(ga_nodeid().eq.0) then
      write(*,*) "got task",i,util_sgroup_mygroup()

      write(*,*) "starting energy calculation",
     >            i,util_sgroup_mygroup()
      end if
      call smd_monomer_energy(i,rtdb)
      if(ga_nodeid().eq.0) then
      write(*,*) "finished energy calculation",
     >            i,util_sgroup_mygroup()
      end if
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
       call util_sgend(rtdb)
       energy = 0.0d0
       if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
      call smd_end(rtdb)
      task_smd_fragment_energy = .true.
      return
      END
