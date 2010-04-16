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
      call util_message("finished smd_geom_monomer_end")

      write(key,'(a,'':converged'')') theory(1:inp_strlen(theory))
      if (.not. rtdb_put(rtdb, key, mt_log, 1, .false.))
     &    call errquit(pname//'  rtdb_put', 0,
     &                 RTDB_ERR)
      call smd_fragment_charge_print()

      smd_monomer_energy = .true.

      call smd_geom_monomer_end()
      call util_message("finished "//pname)
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
      logical oscratch

      pname = "smd_monomer_init"
c
      call smd_message_post("in "//pname)
c
      oscratch = .true.
      call smd_rtdb_push_open_clone(im,oscratch,rtdb)
      
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
      call smd_message_post("out "//pname)

      return
      end subroutine


