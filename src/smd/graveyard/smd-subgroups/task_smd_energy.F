      function task_smd_energy(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include <mpif.h>
#include "global.fh"
      logical task_smd_energy
      character*72 namespace
      character*32 pname
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision eetime4,eetime5,eetime6
      double precision eetime1,eetime2,eetime3
      double precision etime31
      double precision energy,ewald1,ewald2,evdw
      logical operiodic
      integer i,nr,igr
      integer gr_list(1)
      integer gr_size
      character(256) rtdb_fname
      character(256) rtdb_fname_group(2)
      logical master
      integer ierr
c
      integer world_me
c
      pname = "task_smd_energy"
      if(ga_pgroup_get_default().ne.ga_pgroup_get_world())
     $  call errquit('Cannot run SMD in sub-groups yet',0,0)
c
      world_me = ga_nodeid()
c
      call ga_pgroup_sync(ga_pgroup_get_world())
c
c     start the clock
c     ---------------
      etime1=0.0
      etime2=0.0
      etime3=0.0
      etime4=0.0
      etime5=0.0
      etime6=0.0
c
      etime1 = mpi_wtime()
      eetime1 = util_cpusec()
      call smd_init("energy",rtdb)
      etime2 = mpi_wtime()
c
c      write (*,*) 'calling ga_set_all', ga_nodeid()
c      call util_flush(6)
c
      if(ga_nodeid().eq.0) 
     *   call smd_charge_ga_set_all()
c
c      write (*,*) 'done ga_set_all', ga_nodeid()
c      call util_flush(6)
      call ga_sync()

c      write (*,*) 'calling rtdb_close', ga_nodeid()
c      call util_flush(6)

      if(.not.rtdb_close(rtdb, 'keep'))
     +  call errquit('Failed to close group rtdb',0, GEOM_ERR)
c      write (*,*) 'done rtdb_close', ga_nodeid()
c      call util_flush(6)
      call ga_sync()
      call mpi_barrier(mpi_comm_world,ierr)
c
c      write(*,*) "1"//pname,ga_nodeid(), world_me 
c      call util_flush(6)
c      call cpu_time(etime3)
      etime3 = mpi_wtime()
      if(world_me.eq.0) then
c         write(*,*) "1.s"//pname,ga_nodeid(), world_me
c         call util_flush(6)
         call sched_grp_server_code()
      else
c         write(*,*) "1.c"//pname,ga_nodeid(), world_me
c         call util_flush(6)
         call sched_grp_client_code()
      endif
c      write(*,*) "2"//pname,ga_nodeid() 
c      call util_flush(6)
c
c      call cpu_time(etime31)
      etime31 = mpi_wtime()
      call ga_sync()
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c      call cpu_time(etime4)
      etime4 = mpi_wtime()
      call smd_rtdb_get_name(rtdb_fname)
      if (.not. rtdb_open(rtdb_fname,'old',rtdb)) call errquit
     *     (' rtdb_open old failed ',0,0)
c

      etime5 = mpi_wtime()
      eetime5 = util_cpusec()
c
      call smd_end(rtdb)
c
      etime6 = mpi_wtime()
      eetime6 = util_cpusec()
      if(ga_nodeid().eq.0) then
      write(luout,'(/i4,a,2f20.3)') ga_nodeid(),'Total CPU time : ',
     *   (etime6-etime1),(eetime6-eetime1)
      write(luout,'(/i4,a,f20.3)') ga_nodeid(),'Initialization : ',
     *   (etime2-etime1)
      write(luout,'(/i4,a,f20.3)') ga_nodeid(),'Client-Server : ',
     *   (etime31-etime3)
      write(luout,'(/i4,a,f20.3)') ga_nodeid(),'Client-Server+sync : ',
     *   (etime4-etime3)
      write(luout,'(/i4,a,f20.3)') ga_nodeid(),'SMD end : ',
     *   (etime6-etime5)
      end if
c      write(luout,'(/,a,f20.3)')'Verlet         : ',(etime3-etime2)
c      write(luout,'(/,a,f20.3)')'Force          : ',(etime4-etime3)
c
      energy = 0.0d0
      if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
c
      task_smd_energy = .true.
c
      return
      END


      function task_smd_energy2(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
      logical task_smd_energy2
      character*72 namespace
      character*32 pname
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision energy,ewald1,ewald2,evdw
      logical operiodic
      integer i,nr,igr
      integer gr_list(1)
      integer gr_size
      character(256) rtdb_fname
      character(256) rtdb_fname_group(2)
      logical master
c
      logical util_sgstart
      external util_sgstart
      integer util_sgroup_mygroup
      external util_sgroup_mygroup
c
       if(.not.util_sgstart(rtdb))
     $  call errquit('Error in SUBGROUPS init', 0, RTDB_ERR)

      call smd_rtdb_set_handle(rtdb)

      pname = "task_smd_energy"
c      write(*,*) "in "//pname,util_sgroup_mygroup(),rtdb
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
      call smd_init("energy",rtdb)
      call cpu_time(etime2)
c
c      call smd_coords_print(54)
      call smd_vlist_update(.true.,.true.)
      call cpu_time(etime3)
c
      call smd_system_nres(nr)
c
      write(*,*) "number of monomers",nr
      
c      call smd_charge_print(89)
      if (util_sgroup_mygroup().eq.1) then
      call smd_monomer(187)
      end if
c      if (util_sgroup_mygroup().eq.2) then
c      call smd_monomer(74)
c      end if
c      call smd_charge_print(90)
      call cpu_time(etime4)
c
c      write(6,*) "before sync",util_sgroup_mygroup()
c      call util_flush(6)
      call ga_pgroup_sync(ga_pgroup_get_world())
      call smd_end(rtdb)
c

      call cpu_time(etime6)
      write(luout,'(/,a,f20.3)')'Total CPU time : ',(etime6-etime1)
      write(luout,'(/,a,f20.3)')'Initialization : ',(etime2-etime1)
      write(luout,'(/,a,f20.3)')'Verlet         : ',(etime3-etime2)
      write(luout,'(/,a,f20.3)')'Force          : ',(etime4-etime3)

      call util_sgend(rtdb)
      energy = 0.0d0
      if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
      write(*,*) "Goodbye from task_smd",util_sgroup_mygroup()

      task_smd_energy2 = .true.

      return
      END

      function task_smd_energy1(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
      logical task_smd_energy1
      character*72 namespace
      character*32 pname
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision energy,ewald1,ewald2,evdw
      logical operiodic
      integer i,nr,igr
      integer gr_list(1)
      integer gr_size
      character(256) rtdb_fname
      character(256) rtdb_fname_group(2)
      logical master
c
      logical util_sgstart
      external util_sgstart
      integer util_sgroup_mygroup
      external util_sgroup_mygroup
c
       if(.not.util_sgstart(rtdb))
     $  call errquit('Error in SUBGROUPS init', 0, RTDB_ERR)

c      master = ga_nodeid().eq.0
c      call smd_rtdb_get_handle(rtdb)
c      if (.not. rtdb_getfname(rtdb, rtdb_fname)) call
c     >     errquit('rtdb_fname call failed',0,0)
c      rtdb_fname_group(1)="./perm/test001.db"
c      rtdb_fname_group(2)="./perm/test002.db"
c
c      if(ga_nodeid().eq.0) then
c        call util_file_copy(rtdb_fname,rtdb_fname_group(1))
c        call util_file_copy(rtdb_fname,rtdb_fname_group(2))
c
c      end if
c      call ga_sync()
c      write(*,*) "rtdb_fname",rtdb_fname
c
c      if (.not. rtdb_close(rtdb, 'keep')) call errquit(
c     >     ' failed to close rtdb for sg ',0,0)
c
c      gr_list(1) = ga_nodeid()
c      gr_size = 1
c      call ga_sync()
c      igr=ga_pgroup_create(gr_list, gr_size)
c      call ga_pgroup_set_default(igr)
c
c
c      i = ga_nodeid()+1
c      if (.not. rtdb_open(rtdb_fname_group(i),'old',rtdb)) 
c     >     call errquit
c     *     (' rtdb_open old failed ',0,0)
c
c      call ga_sync()
cc      if(util_sgroup_mygroup().eq.1) then
c      if(.true.) then
      call smd_rtdb_set_handle(rtdb)

      pname = "task_smd_energy"
      write(*,*) "in "//pname,util_sgroup_mygroup(),rtdb
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
      call smd_init("energy",rtdb)
      call cpu_time(etime2)
c
c      call smd_coords_print(54)
      call smd_vlist_update(.true.,.true.)
      call cpu_time(etime3)
c
      call smd_system_nres(nr)
c
      write(*,*) "number of monomers",nr
      
c      call smd_charge_print(89)
c      do i=1,1
      if(util_sgroup_mygroup().eq.1) then
      call smd_monomer(1)
      end if
c      end do
c      call smd_charge_print(90)
      call cpu_time(etime4)
c
      call smd_end(rtdb)
c

      call cpu_time(etime6)
      write(luout,'(/,a,f20.3)')'Total CPU time : ',(etime6-etime1)
      write(luout,'(/,a,f20.3)')'Initialization : ',(etime2-etime1)
      write(luout,'(/,a,f20.3)')'Verlet         : ',(etime3-etime2)
      write(luout,'(/,a,f20.3)')'Force          : ',(etime4-etime3)

c      end if
c      if(.not.ga_pgroup_destroy(igr))
c     $     call errquit(pname//'failed to destroy groups', 0, 0)

      call util_sgend(rtdb)
      energy = 0.0d0
      if (.not. rtdb_put(rtdb,'smd:energy',mt_dbl,1,energy))
     $     call errquit(pname//'failed put energy', 0, RTDB_ERR)
      write(*,*) "Goodbye from task_smd",util_sgroup_mygroup()

      task_smd_energy1 = .true.

c      call ga_sync()
      return
      END

      function task_smd_energy0(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
      logical task_smd_energy0
      character*72 namespace
      character*32 pname
      double precision etime1,etime2,etime3,etime4,etime5,etime6
      double precision energy,ewald1,ewald2,evdw
      logical operiodic

      pname = "task_smd_energy"
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
      call smd_init("energy",rtdb)
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
      task_smd_energy0 = .true.

      return
      END
c $Id$
