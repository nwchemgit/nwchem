      subroutine smd_force_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_force,sp_atom
      character*32 tag,pname
      logical result

      pname = "smd_force_init_system"
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "force"
      call smd_system_get_component(sp_force,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      call smd_force_init(sp_force,sp_atom)
c
      return
      end

      subroutine smd_force_init(sp_force,sp_atom)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_force
      character*(*) sp_atom
c
      character*32 pname
      integer na
c
      pname = "smd_force_init"
c
c      write(*,*) "in "//pname
c
c     get total number of atoms 
c     -------------------------
      call smd_atom_ntot(sp_atom,na)
      if(na.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
c     create force data structures
c     ---------------------------
      call smd_namespace_create(sp_force)
      call smd_data_create(sp_force,"force",3*na,MT_DBL)
      return
      end

      subroutine smd_force_print(un)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer un
c
      character*72 sp_force
      character*32 pname,tag
      integer na
      logical result
      integer i_f
      integer i,k
c
      pname = "smd_force_print"
c
      write(*,*) "in "//pname
c
      call smd_system_get_component(sp_force,"force",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no energy ',0,0)

c
      tag = "force"
      call smd_data_get_index(sp_force,tag,i_f,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      call smd_data_get_size(sp_force,tag,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)
      na = na/3

      call smd_util_print_force_array(un,na,dbl_mb(i_f))

      return
      end

      subroutine smd_force_zero()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer un
c
      character*72 sp_force
      character*32 pname,tag
      integer na
      logical result
      integer i_f
      integer i,k
c
      pname = "smd_force_print"
c
      write(*,*) "in "//pname
c
      call smd_system_get_component(sp_force,"force",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no energy ',0,0)

c
      tag = "force"
      call smd_data_get_index(sp_force,tag,i_f,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      call smd_data_get_size(sp_force,tag,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      call dfill(na,
     +             0.d0,
     +             dbl_mb(i_f),
     +             1)


      return
      end

      subroutine smd_force_compute()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_system.fh"
c     
      integer un
c
      character*32 pname
c
      pname = "smd_force_compute"
c
      write(*,*) "in "//pname

      call smd_force_zero()
      call smd_coulomb_compute()
      call smd_pot_compute()
      return
      end
c $Id$
