      logical function task_smd_dynamics(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "smd_system.fh"
      character*72 namespace
      character*32 pname
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision energy,ewald1,ewald2,evdw
      logical operiodic
      integer i
      double precision ke
      integer ndata
      logical vupdate

      pname = "task_smd_dynamics"
      write(*,*) "in "//pname
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
      call smd_init("dynamics",rtdb)
      call cpu_time(etime2)
cc
c      call smd_vel_read(70)
c      call smd_vel_write(71)
cc
c      if (.not.smd_system_ndata(ndata)) 
c     >      call errquit(pname//'failed to get ndata',0,
c     >       RTDB_ERR)
cc
c      write(*,*) "ndata",ndata
c      do i=1,ndata
c        write(93,*) "niter",i
c        call smd_vlist_test(vupdate)
c        call smd_vlist_update(vupdate,.true.)
c        call smd_force_compute()
c        call smd_coords_update()
c        call smd_energy_compute()
c        call smd_energy_print(6)
c      end do
cc
c      call smd_vel_ke_compute(ke)
c      call smd_end(rtdb)
cc
c      if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
c     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
cc
c
c      call cpu_time(etime6)
c      write(luout,'(/,a,f20.3)')'Total CPU time : ',(etime6-etime1)
c      write(luout,'(/,a,f20.3)')'Initialization : ',(etime2-etime1)
c      write(luout,'(/,a,f20.3)')'Verlet         : ',(etime3-etime2)
c      write(luout,'(/,a,f20.3)')'Force          : ',(etime4-etime3)

      write(*,*) "Goodbye from task_smd"
      task_smd_dynamics = .true.

      return
      END
c $Id$
