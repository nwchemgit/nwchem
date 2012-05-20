      function smd_monomer(im)
      implicit none

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "process_task.fh"
#include <mpif.h>
      logical smd_monomer
      integer im
c
      character*72 namespace
      character*32 pname
      integer rtdb
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
      double precision e44,e55
      double precision e1, e2,e3,e4,e5,e6,e7,e8,e9,e10,e11

      double precision util_cpusec
      external util_cpusec
      pname = "smd_monomer"
c      write(*,*) "1 im=",im," in "//pname
      e1 = mpi_wtime()
      call smd_rtdb_get_handle(rtdb)
      if (.not.rtdb_cget(rtdb,"smd:theory",1,theory)) then
       write(0,*) "PLEASE SET smd:theory"
       call errquit('failed to get theory',0,0)
      end if
      call smd_charge_ga_get_all()
c      write(*,*) "1 "//pname
c      call util_flush(6)
      call ga_pgroup_sync(ga_pgroup_get_default())

      call smd_geom_monomer_create(11)
      e2 = mpi_wtime()
      t_gm_crt = t_gm_crt + (e2-e1)
c     write (*,*) "2 im=",im," in "//pname, ga_nodeid()
c     call util_flush(6)

c     write (*,*) "3 im=",im," in "//pname, ga_nodeid()
c     call util_flush(6)
c      call smd_rtdb_get_handle(rtdb)
c      write (*,*) "4 im=",im," in "//pname, ga_nodeid()
c      call util_flush(6)
      e4 = mpi_wtime()
      e44=util_cpusec()
      if(.not.rtdb_put(rtdb,'smd:im',mt_int,1,im))
     +     call errquit('rtdb_put smd:im failed',0,0)
      result = task_energy_doit(rtdb,theory,qm_energy)
      e5 = mpi_wtime()
      e55=util_cpusec()
      if(ga_nodeid().eq.0) write(*,*)"time energy", (e5-e4),(e55-e44)
      t_doit = t_doit + (e5-e4)
c      write (*,*) "5 im=",im," in "//pname, ga_nodeid()
c      call util_flush(6)
      e6 = mpi_wtime()
      if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
     +     call errquit('esp: rtdb_put newgrd failed',0,0)
c      write (*,*) "6 im=",im," in "//pname, ga_nodeid()
c      call util_flush(6)

      if (theory.eq.'dft') then
       tag='dft:output vectors'
      else
       tag='scf:output vectors'
      endif
      if (.not.rtdb_cget(rtdb,tag,1,
     >                   movecs)) then
       write (6,*) 'error in rtdb_cget'
      end if
c
      if (.not.rtdb_cput(rtdb,'esp:input vectors',1,
     >                   movecs)) then
       write (6,*) 'error in rtdb_cput'
      end if
c
      e7 = mpi_wtime()
      t_rtdb = t_rtdb + (e7-e6)
c
c      if(.not.esp(rtdb))
c     &    call errquit(pname//'failed esp call', 100,
c     &                 RTDB_ERR)
      e8 = mpi_wtime()
      t_esp = t_esp + (e8-e7)
c      write (*,*) "7 im=",im," in "//pname, ga_nodeid()
c      call util_flush(6)

c      call smd_geom_push_charges(im)
      e9 = mpi_wtime()
      t_gm_push = t_gm_push + (e9-e8)

c      if(ga_nodeid().eq.0)
c     &     call smd_charge_ga_set1(im)
c     call smd_charge_ga_set(im)
      e10 = mpi_wtime()
      t_chrg_set = t_chrg_set + (e10-e9)
 
      call smd_geom_monomer_end()
      e11 = mpi_wtime()
      t_im_end = t_im_end + (e11-e10)

c      write(*,*) "Goodbye from "//pname,util_sgroup_mygroup(),qm_energy
c      call util_flush(6)
      if(ga_nodeid().eq.0)
     >  call util_message("smd monomer done")
      smd_monomer = .true.

      return
      end function

c $Id$
