      function smd_monomer_energy(im,theory)
      implicit none
c     calculates monomer (im) energy
c     in the field (classical) of others
c     im[IN]  monomer id
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
      logical smd_monomer_energy
      integer im
      character*(*) theory
c
      integer rtdb
      character*72 operation
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
      character*256 rtdb_name
      double precision e1, e2,e3,e4,e5,e6,e7,e8,e9,e10,e11
      character*255 prefix, prefix0,key
      integer i1
      integer group_id

      pname = "smd_monomer_energy"

      call smd_monomer_init(im,theory,operation,rtdb)
      result = task_energy_doit(rtdb,theory,qm_energy)
c
      if(ga_nodeid().eq.0)
     >  write(*,*) "finished energy",theory,qm_energy
c
      if(theory.ne."dft") goto 100
      if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
     +     call errquit('esp: rtdb_put newgrd failed',0,0)


      if(.not.esp(rtdb))
     &    call errquit(pname//'failed esp call', 100,
     &                 RTDB_ERR)
c
      call util_message("finished esp")
c
      call smd_geom_push_charges(im)
      if(ga_nodeid().eq.0)
     &     call smd_fragment_charge_set(im)

      call util_message("finished smd_geom_push_charges")

100   continue
      call smd_geom_monomer_end()
      call util_message("smd_geom_monomer_end")

      write(key,'(a,'':converged'')') theory(1:inp_strlen(theory))
      if (.not. rtdb_put(rtdb, key, mt_log, 1, .false.))
     &    call errquit(pname//'  rtdb_put', 0,
     &                 RTDB_ERR)
      call smd_fragment_charge_print()

      call smd_rtdb_close("delete")
      call util_file_prefix_set(prefix0)
      call smd_rtdb_open_master(rtdb)
      smd_monomer_energy = .true.

      return
      end function
      function smd_monomer_energy2(im,theory)
      implicit none
c     calculates monomer (im) energy
c     in the field (classical) of others
c     im[IN]  monomer id
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
      logical smd_monomer_energy2
      integer im
      character*(*) theory
c
      integer rtdb
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
      character*256 rtdb_name
      double precision e1, e2,e3,e4,e5,e6,e7,e8,e9,e10,e11
      character*255 prefix, prefix0,key
      integer i1
      integer group_id

      pname = "smd_monomer_energy"
c
c
c     set prefix for all generated files (e.g. movecs)
c     so that it does not clash with other calculations
c     -------------------------------------------------
      call util_file_prefix_get(prefix0)
c      i1 = inp_strlen(prefix0)
c      write(prefix,'(A,I4.4)')
c     >    prefix0(1:i1),im
c      call util_file_prefix_set(prefix)
c      call util_message("prefix is set to "//prefix)
c     -------------------------------------------------
c     get rtdb handle from smd database (in-core)
c      call smd_rtdb_get_handle(rtdb)
      call smd_global_update("coords","get")
      call smd_rtdb_close("keep")
c      call smd_rtdb_open_clone(im,rtdb)
      
c      call smd_rtdb_get_string("smd:rtdb_name",1,
c     >                           rtdb_name,result)

c     -------------------------------------------------
c     set the theory (but this may not be necessary)
c     if I do task fragment dft energy
      if (.not.(rtdb_cput(rtdb,'task:theory',1,theory)))
     *     call errquit('rtdb_put task:theory failed',0,0)
c
c     pull the global charges onto local set
c      should this be done outside this sub
      call smd_fragment_charge_get_all()
c
c      create geometry object
c      ---------------------
      call smd_geom_monomer_create(im)
      if(ga_nodeid().eq.0)
     >  write(*,*) "smd theory is set to",theory,group_id
      result = task_energy_doit(rtdb,theory,qm_energy)
c
      if(ga_nodeid().eq.0)
     >  write(*,*) "finished energy",theory,qm_energy
c
      if(theory.ne."dft") goto 100
      if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
     +     call errquit('esp: rtdb_put newgrd failed',0,0)


      if(.not.esp(rtdb))
     &    call errquit(pname//'failed esp call', 100,
     &                 RTDB_ERR)
c
      call util_message("finished esp")
c
      call smd_geom_push_charges(im)
      if(ga_nodeid().eq.0)
     &     call smd_fragment_charge_set(im)

      call util_message("finished smd_geom_push_charges")

100   continue
      call smd_geom_monomer_end()
      call util_message("smd_geom_monomer_end")

      write(key,'(a,'':converged'')') theory(1:inp_strlen(theory))
      if (.not. rtdb_put(rtdb, key, mt_log, 1, .false.))
     &    call errquit(pname//'  rtdb_put', 0,
     &                 RTDB_ERR)
      call smd_fragment_charge_print()

      call smd_rtdb_close("delete")
      call util_file_prefix_set(prefix0)
      call smd_rtdb_open_master(rtdb)
      smd_monomer_energy2 = .true.

      return
      end function

      subroutine smd_monomer_init(im,theory,operation,rtdb)
      implicit none
c     calculates monomer (im) energy
c     in the field (classical) of others
c     im[IN]  monomer id
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
      integer im
      character*(*) theory
      character*(*) operation
      integer rtdb
c
      character*72 namespace
      character*32 pname
      integer group_id

      pname = "smd_monomer_init"
c
      call smd_rtdb_close("keep")
c      call smd_rtdb_open_clone(im,rtdb)
      
c     set the theory 
c     --------------
      if (.not.(rtdb_cput(rtdb,'task:theory',1,theory)))
     *     call errquit('rtdb_put task:theory failed',0,0)
c
c     pull the global charges/coords onto local set
      call smd_fragment_charge_get_all()
      call smd_global_update("coords","get")
c
c      create geometry object
c      ---------------------
      call smd_geom_monomer_create(im)
c
      call util_message("finished "//pname)

      return
      end subroutine

c      subroutine smd_monomer_end(im,theory,operation,rtdb)
c      implicit none
c#include "inp.fh"
c#include "mafdecls.fh"
c#include "rtdb.fh"
c#include "stdio.fh"
c#include "errquit.fh"
c#include "global.fh"
c      integer im
c      character*(*) theory
c      character*(*) operation
cc
c      integer rtdb
c      character*72 namespace
c      character*32 pname
c      logical result
c      double precision qm_energy
c      logical mode
c      character*72 tag
c      character*256 movecs
c      character*256 rtdb_name
c      double precision e1, e2,e3,e4,e5,e6,e7,e8,e9,e10,e11
c      character*255 prefix, prefix0,key
c      integer i1
c      integer group_id
c
c      pname = "smd_monomer_end"
c
c      call smd_geom_monomer_end()
c
c      call smd_rtdb_close("delete")
c      call util_file_prefix_set(prefix0)
c      call smd_rtdb_open_master(rtdb)
c      smd_monomer_energy = .true.
c
c      return
c      end subroutine

      function smd_monomer_energy1(im)
      implicit none
c     calculates monomer (im) energy
c     in the field (classical) of others
c     im[IN]  monomer id
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
      logical smd_monomer_energy1
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
      character*255 prefix, prefix0,key
      integer i1
      integer group_id

      pname = "smd_monomer_energy1"
c
c
      call util_file_prefix_get(prefix0)
      i1 = inp_strlen(prefix0)
      write(prefix,'(A,I4.4)')
     >    prefix0(1:i1),im
      call util_file_prefix_set(prefix)
c
      if (.not.(rtdb_cget(rtdb,'smd:theory',1,theory)))
     *     theory = "dft"
      if (.not.(rtdb_cput(rtdb,'task:theory',1,theory)))
     *     call errquit('rtdb_put task:theory failed',0,0)

      call smd_fragment_charge_get_all()
c
      call smd_geom_monomer_create(im)
      if(ga_nodeid().eq.0)
     >  write(*,*) "smd theory is set to",theory,group_id
      result = task_energy_doit(rtdb,theory,qm_energy)
c
      if(ga_nodeid().eq.0)
     >  write(*,*) "finished energy",theory,group_id
c
      if(theory.ne."dft") goto 100
      if(.not.rtdb_put(rtdb,'esp:newgrd',mt_int,1,1))
     +     call errquit('esp: rtdb_put newgrd failed',0,0)


      if(.not.esp(rtdb))
     &    call errquit(pname//'failed esp call', 100,
     &                 RTDB_ERR)
c
      call util_message("finished esp")
c
      call smd_geom_push_charges(im)
      if(ga_nodeid().eq.0)
     &     call smd_fragment_charge_set(im)

      call util_message("finished smd_geom_push_charges")

100   continue
      call smd_geom_monomer_end()
      call util_message("smd_geom_monomer_end")

      write(key,'(a,'':converged'')') theory(1:inp_strlen(theory))
      if (.not. rtdb_put(rtdb, key, mt_log, 1, .false.))
     &    call errquit(pname//'  rtdb_put', 0,
     &                 RTDB_ERR)
      call smd_fragment_charge_print()

      call util_file_prefix_set(prefix0)
      smd_monomer_energy1 = .true.

      return
      end function

