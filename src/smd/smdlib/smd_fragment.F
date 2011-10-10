      subroutine smd_fragment_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_fragment,sp_param,sp_atom,sp_type
      character*32 tag,pname
      logical result

      pname = "smd_charge_init_system"
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "parameters"
      call smd_system_get_component(sp_param,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "fragment"
      call smd_system_get_component(sp_fragment,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "type"
      call smd_system_get_component(sp_type,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      call ga_sync()
      call smd_fragment_init(sp_fragment)
c      call smd_fragment_charge_print()
      call ga_sync()
c
      return
      end

      subroutine smd_fragment_init(sp_fragment)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_charge_data.fh"
c     
      character*(*) sp_fragment
c
      character*32 pname
      integer na
c
      integer ndims, type, ga_charges
      integer ga_energies,ga_charges0
      integer dims(1), chunk(1)
      integer rtdb
c
      pname = "smd_charge_init"
c
c      write(*,*) "in "//pname
c
      type = MT_DBL
      ndims = 1
c
c     get total number of atoms 
c     -------------------------
      call smd_atom_ntot(na)
      if(na.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
c     create data structures
c     ---------------------------
      call smd_namespace_create(sp_fragment)
      if(.not.nga_create(mt_dbl, 1, na, 
     >    'fragment_charges',1, ga_charges)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)
      call smd_data_put_int(sp_fragment,"charge",1,ga_charges)
c
      if(.not.ga_duplicate(ga_charges,ga_charges0,"fragment_charges0"))
     >  call errquit(
     >       pname//'ga_duplicate error',0, RTDB_ERR)
      call smd_data_put_int(sp_fragment,"charge0",1,ga_charges0)
c
      if(.not.nga_create(mt_dbl, 1, na, 
     >    'fragment_energies',1, ga_energies)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)
      call smd_data_put_int(sp_fragment,"energy",1,ga_energies)

      call ga_sync()
      if(ga_nodeid().eq.0)
     > call smd_fragment_charge_set_all()
      call ga_sync()
     
      call ga_copy(ga_charges,ga_charges0)

      return
      end

      subroutine smd_fragment_charge_set_all1()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer nr
      integer i
      character*36 pname
      logical mode
c
      pname = "smd_charge_set_all"

      call smd_system_nres(nr)

      do i=1,nr
       call smd_fragment_charge_set(i)
      end do

      return
      end

      subroutine smd_fragment_charge_set_all()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 pname,tag
      logical result
      integer i_q
      integer ga_charges
      integer na
      integer lo(1), hi(1), ld(1)
c
      pname = "smd_fragment_charge_set_all"
c
      tag = "charge"
      call smd_get_ind_dim(tag,i_q,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      call smd_fragment_charge_handle(ga_charges)

      lo(1) = 1
      hi(1) = na
      ld(1) = 1
      call nga_put(ga_charges,lo, hi, dbl_mb(i_q), ld)
      return
      end

      subroutine smd_fragment_charge_get_all()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 pname,tag
      logical result
      integer i_q
      integer ga_charges
      integer na
      integer lo(1), hi(1), ld(1)
c
      pname = "smd_fragment_charge_get_all"
c
      tag = "charge"
      call smd_get_ind_dim(tag,i_q,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      call smd_fragment_charge_handle(ga_charges)

      lo(1) = 1
      hi(1) = na
      ld(1) = 1
      call nga_get(ga_charges,lo, hi, dbl_mb(i_q), ld)
      return
      end

      subroutine smd_fragment_charge_set(im)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer im
c
      character*32 pname,tag
      logical result
      integer i_q,i_q0,i_id
      integer i,j
      double precision q0
      integer rtdb,ga_charges
      integer na,nm,i0
      integer lo(1), hi(1), ld(1)
      integer i_ir
c
      pname = "smd_fragment_charge_set"
c
      tag = "charge"
      call smd_get_ind_dim(tag,i_q,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      tag = "atom:resid"
      call smd_get_ind(tag,i_ir,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      call smd_fragment_charge_handle(ga_charges)

      call smd_geom_monomer_index_size(im,na,int_mb(i_ir),nm,i0)

      lo(1) = i0
      hi(1) = i0+nm-1
      ld(1) = 1
      call nga_put(ga_charges,lo, hi, dbl_mb(i_q+i0-1), ld)
      return
      end
      subroutine smd_fragment_charge_get_all1()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer nr
      integer i
      character*32 pname
c
      pname = "smd_charge_get_all"
c
c      write(*,*) "in "//pname

      call smd_system_nres(nr)

      do i=1,nr
       call smd_fragment_charge_get(i)
      end do

      return
      end

      subroutine smd_fragment_charge_get(im)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer im
      character*32 sp_fragment
      character*32 sp_atom
c
      character*32 pname,tag
      logical result
      integer i_q,i_q0,i_id
      integer i,j
      double precision q0
      integer rtdb,ga_charges
      integer na,nm,i0
      integer lo(1), hi(1), ld(1)
      integer i_ir
c
      pname = "smd_charge_get"
c
c      write(*,*) "in "//pname
c
c
      tag = "charge"
      call smd_get_ind_dim(tag,i_q,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      tag = "atom:resid"
      call smd_get_ind(tag,i_ir,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
 
      call smd_geom_monomer_index_size(im,na,int_mb(i_ir),nm,i0)
      call smd_fragment_charge_handle(ga_charges)

      lo(1) = i0
      hi(1) = i0+nm-1
      ld(1) = 1
      call nga_get(ga_charges,lo, hi, dbl_mb(i_q+i0-1), ld)

      return
      end

      subroutine smd_fragment_charge_print()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 pname,tag
      logical result
      integer rtdb,ga
c
      pname = "smd_charge_print"
c

      call smd_fragment_charge_handle(ga)
c      call ga_print_distribution(ga)
c      call ga_print(ga)

c      call smd_fragment_energy_handle(ga)
c      call ga_print_distribution(ga)
c      call ga_print(ga)

      return
      end

      subroutine smd_fragment_charge_handle(ga_handle)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "smd_const_data.fh"
c     
      integer ga_handle
c     
      logical result
      character*(smd_string_size) tag,sp,pname
c
      pname  = "smd_fragment_charge_handle"
c
      tag = "fragment"
      call smd_system_get_component(sp,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_data_get_int(sp,"charge",1,ga_handle,result)
      if(.not.result)
     >  call errquit(
     >       pname//' getting handle'//tag,0,0)
      return
      end

      subroutine smd_fragment_charge_diff(eps)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "smd_const_data.fh"
      double precision eps
c     
      integer ga
      integer ga0
      integer p_world,p_default
c     
      logical result
      character*(smd_string_size) tag,sp,pname
c
      pname  = "smd_fragment_charge_handle"
c
      p_world   = ga_pgroup_get_world()
      p_default = ga_pgroup_get_default()
      call ga_pgroup_sync(p_world)
c
      call ga_pgroup_set_default(p_world)
c
      tag = "fragment"
      call smd_system_get_component(sp,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_data_get_int(sp,"charge",1,ga,result)
      if(.not.result)
     >  call errquit(
     >       pname//' getting handle'//tag,0,0)
c
      call smd_data_get_int(sp,"charge0",1,ga0,result)
      if(.not.result)
     >  call errquit(
     >       pname//' getting handle'//tag,0,0)

      call  ga_add(1.0d0,ga,-1.0d0,ga0,ga0)
      call ga_norm_infinity(ga0,eps) 
      if(ga_nodeid().eq.0) then
        write(*,*) "current charge difference",eps
        call util_flush(6)
      end if

      call ga_copy(ga,ga0)
     
      call ga_pgroup_set_default(p_default)

      return
      end

      subroutine smd_fragment_energy_handle(ga_handle)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "smd_const_data.fh"
c     
      integer ga_handle
c     
      logical result
      character*(smd_string_size) tag,sp,pname
c
      pname  = "smd_fragment_charge_handle"
c
      tag = "fragment"
      call smd_system_get_component(sp,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_data_get_int(sp,"energy",1,ga_handle,result)
      if(.not.result)
     >  call errquit(
     >       pname//' getting handle'//tag,0,0)
      return
      end
c $Id$
