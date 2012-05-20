      function task_smd_energy(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
      logical task_smd_energy
      logical ogroups,result,ofragment
c
      logical task_smd_energy_groups
      external task_smd_energy_groups
      logical task_smd_fragment
      external task_smd_fragment
      logical task_smd_energy_sp
      external task_smd_energy_sp
c
      double precision t1,t2
c
      t1 = util_cpusec()
      call smd_group_inquire(rtdb,ogroups)
      if (.not.rtdb_get(rtdb, 'task:fragment', mt_log, 1,ofragment))
     +   ofragment = .false.
c
      if(ofragment) then
       result = task_smd_fragment(rtdb)
       goto 100
      end if
c
      call smd_init("energy",rtdb)
c
      if(ogroups) then
c        result = task_smd_energy_groups(rtdb)
      else
        result = task_smd_energy_sp(rtdb)
      end if

100   continue

      task_smd_energy = .true.

      t2 = util_cpusec()
      if(ga_nodeid().eq.0) write(*,*) "total time",t2-t1
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
c
c      function task_smd_energy_groups(rtdb)
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
c      logical task_smd_energy_groups
c      character*72 namespace
c      character*32 pname
c      double precision etime1,etime2,etime3,etime4,etime5,etime6
c      double precision energy,ewald1,ewald2,evdw
c      logical operiodic
c      integer i
cc
c      integer util_sgroup_mygroup
c      external util_sgroup_mygroup
cc
c      logical smd_tasklist_get
c      external smd_tasklist_get
cc
c      integer util_sgroup_numgroups
c      external util_sgroup_numgroups
c
c      pname = "task_smd_energy_groups"
cc
c      call smd_tasklist_init()
cc     
c      call ga_sync()
cc
c      call smd_init_group(rtdb)
cc
cc      write(*,*) "number of groups",util_sgroup_numgroups(),
cc:     >           util_sgroup_mygroup()
c
c99    continue
c      if(.not.smd_tasklist_get(i))
c     >  goto 100 
c      if(ga_nodeid().eq.0) then
c      write(*,*) "got task",i,util_sgroup_mygroup()
c
c      write(*,*) "starting energy calculation",
c     >            i,util_sgroup_mygroup()
c      end if
c      call smd_monomer_energy(i)
c      if(ga_nodeid().eq.0) then
c      write(*,*) "finished energy calculation",
c     >            i,util_sgroup_mygroup()
c      end if
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
c
c      call ga_pgroup_sync(ga_pgroup_get_world())
c       call util_sgend(rtdb)
c       energy = 0.0d0
c       if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
c     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
c      call smd_end(rtdb)
c      task_smd_energy_groups = .true.
c      return
c      END
c $Id$
