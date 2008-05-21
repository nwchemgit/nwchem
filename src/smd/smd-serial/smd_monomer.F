      function smd_monomer_energy(im)
      implicit none

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
      logical smd_monomer_energy
      integer im,rtdb
c
      character*72 namespace
      character*32 pname
      logical result
      double precision qm_energy
      logical task_energy_doit
      external task_energy_doit
      logical esp
      external esp
      integer util_sgroup_mygroup
      external util_sgroup_mygroup
      logical mode
      character*72 tag
      character*256 movecs
      character*16 theory
      double precision e1, e2,e3,e4,e5,e6,e7,e8,e9,e10,e11

      pname = "smd_monomer"
      theory = "dft"
c
c
      call smd_group_rtdb_clone(rtdb) 
c
      if (.not.(rtdb_cput(rtdb,'task:theory',1,theory)))
     *     call errquit('rtdb_put task:theory failed',0,0)

c
      call smd_geom_monomer_create(im)
      result = task_energy_doit(rtdb,theory,qm_energy)
c
      if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
     +     call errquit('esp: rtdb_put newgrd failed',0,0)

c      if (theory.eq.'dft') then
c       tag='dft:output vectors'
c      else
c       tag='scf:output vectors'
c      endif
c      if (.not.rtdb_cget(rtdb,tag,1,
c     >                   movecs)) then
c       write (6,*) 'error in rtdb_cget'
c      end if
c
c
c      if (.not.rtdb_cput(rtdb,'esp:input vectors',1,
c     >                   movecs)) then
c       write (6,*) 'error in rtdb_cput'
c      end if
c
      if(.not.esp(rtdb))
     &    call errquit(pname//'failed esp call', 100,
     &                 RTDB_ERR)
c
c      call smd_geom_push_charges(im)

      call smd_geom_monomer_end()

      call smd_group_rtdb_release(rtdb)
      smd_monomer_energy = .true.

      return
      end function
c      function smd_monomer(im)
c      implicit none
c
c#include "inp.fh"
c#include "mafdecls.fh"
c#include "rtdb.fh"
c#include "stdio.fh"
c#include "errquit.fh"
c#include "global.fh"
c#include "process_task.fh"
c#include <mpif.h>
c      logical smd_monomer
c      integer im
cc
c      character*72 namespace
c      character*32 pname
c      integer rtdb
c      logical result
c      double precision qm_energy
c      logical task_energy_doit
c      external task_energy_doit
c      logical esp
c      external esp
c      integer util_sgroup_mygroup
c      external util_sgroup_mygroup
c      logical mode
c      character*72 tag
c      character*256 movecs
c      character*16 theory
c      double precision e1, e2,e3,e4,e5,e6,e7,e8,e9,e10,e11
c
c      pname = "smd_monomer"
cc      write(*,*) "1 im=",im," in "//pname
c      e1 = mpi_wtime()
c      call smd_rtdb_get_handle(rtdb)
c      if (.not.rtdb_cget(rtdb,"smd:theory",1,theory)) then
c       write(0,*) "PLEASE SET smd:theory"
c       call errquit('failed to get theory',0,0)
c      end if
cc      call smd_charge_ga_get_all()
cc      write(*,*) "1 "//pname
cc      call util_flush(6)
c      call ga_pgroup_sync(ga_pgroup_get_default())
c
c      call smd_geom_monomer_create(im)
c      e2 = mpi_wtime()
c      t_gm_crt = t_gm_crt + (e2-e1)
cc     write (*,*) "2 im=",im," in "//pname, ga_nodeid()
cc     call util_flush(6)
c
cc     write (*,*) "3 im=",im," in "//pname, ga_nodeid()
cc     call util_flush(6)
cc      call smd_rtdb_get_handle(rtdb)
cc      write (*,*) "4 im=",im," in "//pname, ga_nodeid()
cc      call util_flush(6)
c      e4 = mpi_wtime()
c      result = task_energy_doit(rtdb,theory,qm_energy)
c      e5 = mpi_wtime()
c      t_doit = t_doit + (e5-e4)
c      if(ga_nodeid().eq.0) then
c        write(*,*) "time energy monomer",im, e5-e4
c        write(*,*) "time geom monomer",im,e2-e1
c      end if
c      call util_message("task_energy_doit done")
c      e6 = mpi_wtime()
c      if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
c     +     call errquit('esp: rtdb_put newgrd failed',0,0)
cc      write (*,*) "6 im=",im," in "//pname, ga_nodeid()
cc      call util_flush(6)
c
c      if (theory.eq.'dft') then
c       tag='dft:output vectors'
c      else
c       tag='scf:output vectors'
c      endif
c      if (.not.rtdb_cget(rtdb,tag,1,
c     >                   movecs)) then
c       write (6,*) 'error in rtdb_cget'
c      end if
cc
c      if (.not.rtdb_cput(rtdb,'esp:input vectors',1,
c     >                   movecs)) then
c       write (6,*) 'error in rtdb_cput'
c      end if
cc
c      e7 = mpi_wtime()
c      t_rtdb = t_rtdb + (e7-e6)
cc
cc      if(.not.esp(rtdb))
cc     &    call errquit(pname//'failed esp call', 100,
cc     &                 RTDB_ERR)
c      e8 = mpi_wtime()
c      t_esp = t_esp + (e8-e7)
cc      write (*,*) "7 im=",im," in "//pname, ga_nodeid()
cc      call util_flush(6)
c
cc      call smd_geom_push_charges(im)
c      e9 = mpi_wtime()
c      t_gm_push = t_gm_push + (e9-e8)
c
c      call util_message("before call smd_charge_ga_set1")
cc      if(ga_nodeid().eq.0)
cc     &     call smd_charge_ga_set1(im)
cc     call smd_charge_ga_set(im)
c      call util_message("after call smd_charge_ga_set1")
c      e10 = mpi_wtime()
c      t_chrg_set = t_chrg_set + (e10-e9)
c 
c      call smd_geom_monomer_end()
c      call util_message("after smd_geom_monomer_end")
c      e11 = mpi_wtime()
c      t_im_end = t_im_end + (e11-e10)
c
cc      write(*,*) "Goodbye from "//pname,util_sgroup_mygroup(),qm_energy
cc      call util_flush(6)
c      smd_monomer = .true.
c
c      return
c      end function
c
