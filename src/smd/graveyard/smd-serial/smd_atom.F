      subroutine smd_atom_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_atom
      character*32 tag,pname
      logical result

      pname = "smd_atom_init_system"
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_atom_init(sp_atom)

      return
      end

      subroutine smd_atom_init(namespace)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "geom.fh"
c     
      character*(*) namespace
c
      character*32 pname
      character*80 tag
      integer nt
      integer i_t,i_tr,i_ir,i_qf
      integer i_iconst
c
      logical geom_tag_to_charge_gen
      external geom_tag_to_charge_gen
c
      pname = "smd_atom_init"
c
      write(*,*) "in "//pname
c
c     get total number of atoms from pdb file
c     ---------------------------------------
      call smd_coordfile_natoms(nt)
      if(nt.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
c
c     create atom namespace
c     ---------------------
      call smd_namespace_create(namespace)
c
c     create atom data structures
c     ---------------------------
      tag = "atom:iconst"
      call smd_data_create_get(namespace,tag,1,MT_INT,i_iconst)
      tag = "atom:name"
      call smd_data_create_get(namespace,tag,nt,MT_STR,i_t)
      tag  = "atom:resname"
      call smd_data_create_get(namespace,tag,nt,MT_STR,i_tr)
      tag = "atom:resid"
      call smd_data_create_get(namespace,tag,nt,MT_INT,i_ir)
      tag = "atom:formal_charge"
      call smd_data_create_get(namespace,tag,nt,MT_DBL,i_qf)

c
c    don't ask
c    --------
      int_mb(i_iconst) = nt
c
c     fill in names from pdb file
c     ---------------------------
      call smd_coordfile_read_atomres(nt,
     +                         byte_mb(i_t),
     +                         byte_mb(i_tr),
     +                         int_mb(i_ir))

c
      if(.not.geom_tag_to_charge_gen(nt,smd_string_size,
     >                        byte_mb(i_t),dbl_mb(i_qf)))
     >  call errquit(
     >       pname//'error setting formal charge',0, RTDB_ERR)

      return
      end

      subroutine smd_atom_ntot(ntot)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer ntot
c
      character*72 tag
      character*30 pname
      integer i_iconst
      logical result

      pname = "smd_atom_ntot"
      tag = "atom:iconst"
      call smd_get_ind(tag,i_iconst,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting ntot '//tag,0, RTDB_ERR)
      ntot = int_mb(i_iconst)

      return
      end

      subroutine smd_atom_nres(namespace,nr)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) namespace
      integer nr
c
      character*72 tag
      character*30 pname
      integer i_ir
      logical result
      integer nt
      integer ir,ir0,i

      pname = "smd_atom_nres"
      tag = "atom:resid"
      call smd_data_get_index(namespace,tag,i_ir,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting '//tag,0, RTDB_ERR)
      call smd_data_get_size(namespace,tag,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting size '//tag,0, RTDB_ERR)

      nr = 0
      ir0 = 0
      do i=1,nt
       ir = int_mb(i_ir+i-1)
       if(ir0.ne.ir) then
          ir0=ir
          nr = nr + 1
       end if
      end do
       
      return
      end

      subroutine smd_atom_ntot1(namespace,ntot)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) namespace
      integer ntot
c
      character*72 tag
      character*30 pname
      logical result

      pname = "smd_atom_ntot"
      tag = "atom:xyz"
      call smd_data_get_size(namespace,tag,ntot,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting size for '//tag,0, RTDB_ERR)
      ntot = ntot/3

      return
      end

c $Id$
