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
      logical result,ofragment
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
      if (.not.rtdb_get(rtdb, 'task:fragment', mt_log, 1,ofragment))
     +   ofragment = .false.
c
      if(ofragment) then
       call util_message("executing fragment calculation")
       result = task_smd_fragment(rtdb)
      else 
        result = task_smd_energy_sp(rtdb)
      end if

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
      call smd_init("energy",rtdb)
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

c $Id$
