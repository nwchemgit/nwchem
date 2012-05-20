      function smd_monomer_energy(im)
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

      pname = "smd_monomer_energy"
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
      smd_monomer_energy = .true.

      return
      end function

c $Id$
