      function smd_monomer(im)
      implicit none

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
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


      pname = "smd_monomer"
      write(*,*) "in "//pname
c
      theory = "dft"
      call smd_geom_monomer_create(im)

      call ga_pgroup_sync(ga_pgroup_get_default())

      call smd_rtdb_get_handle(rtdb)
      result = task_energy_doit(rtdb,theory,qm_energy)

        if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
     +     call errquit('esp: rtdb_put newgrd failed',0,0)

     
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
      if(.not.esp(rtdb))
     &    call errquit(pname//'failed esp call', 100,
     &                 RTDB_ERR)

      call smd_geom_push_charges(im)

      if(ga_nodeid().eq.0)
     & call smd_charge_ga_set1(im)
c      call smd_charge_ga_set(im)
 
 
      call smd_geom_monomer_end()

      write(*,*) "Goodbye from "//pname,util_sgroup_mygroup(),qm_energy
      smd_monomer = .true.

      return
      end function
c $Id$
