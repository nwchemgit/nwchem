      function task_smd_energy(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      logical task_smd_energy
      logical ogroups,result
c
      logical task_smd_energy_groups
      external task_smd_energy_groups
      logical task_smd_energy_sp
      external task_smd_energy_sp
c
      call smd_init("energy",rtdb)

      call smd_groups_inquire(rtdb,ogroups)
      if(ogroups) then
        result = task_smd_energy_groups(rtdb)
      else
        result = task_smd_energy_sp(rtdb)
      end if

      write(*,*) "Goodbye from task_smd"
      task_smd_energy = .true.

      return
      END

      function task_smd_energy_sp(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      logical task_smd_energy_sp
      character*72 namespace
      character*32 pname
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision energy,ewald1,ewald2,evdw
      logical operiodic

      pname = "task_smd_energy"
c      write(*,*) "in "//pname
c
c     start the clock
c     ---------------
      etime1=0.0
      etime2=0.0
      etime3=0.0
      etime4=0.0
      etime5=0.0
      etime6=0.0
      call cpu_time(etime1)
c
      call cpu_time(etime1)
      call cpu_time(etime2)
c
      call smd_vlist_update(.true.,.true.)
      call cpu_time(etime3)
c
      call smd_force_compute()
      call cpu_time(etime4)
c
      call smd_energy_compute()
      call cpu_time(etime5)
c
      call smd_energy_print(6)
c
      call smd_end(rtdb)
c
      if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
c

      call cpu_time(etime6)
      write(luout,'(/,a,f20.3)')'Total CPU time : ',(etime6-etime1)
      write(luout,'(/,a,f20.3)')'Initialization : ',(etime2-etime1)
      write(luout,'(/,a,f20.3)')'Verlet         : ',(etime3-etime2)
      write(luout,'(/,a,f20.3)')'Force          : ',(etime4-etime3)

      write(*,*) "Goodbye from task_smd"
      task_smd_energy_sp = .true.

      return
      END

      function task_smd_energy_groups(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
      logical task_smd_energy_groups
      character*72 namespace
      character*32 pname
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision energy,ewald1,ewald2,evdw
      logical operiodic
      integer i
c
      integer util_sgroup_mygroup
      external util_sgroup_mygroup
c
      logical smd_tasklist_get
      external smd_tasklist_get
c
      integer util_sgroup_numgroups
      external util_sgroup_numgroups

      pname = "task_smd_energy_groups"
c      write(*,*) "in "//pname
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
      call smd_monomer_energy(i)
      if(ga_nodeid().eq.0) then
      write(*,*) "finished energy calculation",
     >            i,util_sgroup_mygroup()
      end if
c      call smd_vlist_update(.true.,.true.)
cc
c      call smd_force_compute()
cc
c      call smd_energy_compute()
cc
c      call smd_energy_print(6)
c
c
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
      task_smd_energy_groups = .true.
      return
      END
