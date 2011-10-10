      subroutine smd_param_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_param
      character*32 tag,pname
      logical result

      pname = "smd_param_init_system"
c
      tag = "parameters"
      call smd_system_get_component(sp_param,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_param_init(sp_param)

      return
      end

      subroutine smd_param_init(namespace)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
c     
      character*(*) namespace
c
      character*32 pname
      character*80 tag
      character*255 filename
      integer na
      integer i_ta,i_tt,i_tr,i_q,i_m
      logical result
c
      pname = "smd_param_init"
c
c      write(*,*) "in "//pname
c
c
c     get total number of atoms in parameter file
c     -------------------------------------------
      call smd_paramfile_natoms(na)
c      write(*,*) "number of atoms in parameter space",na
      if(na.le.0)
     >  call errquit(
     >       pname//'no parameters fond',0, 0)

c
      call smd_namespace_create(namespace)
c
c     create parameter data structures
c     ---------------------------
      tag = "param:atomname"
c      call smd_data_create_get(namespace,tag,16*na,MT_BYTE,i_ta)
      call smd_data_create_get(namespace,tag,na,MT_STR,i_ta)
      tag = "param:typename"
c      call smd_data_create_get(namespace,tag,16*na,MT_BYTE,i_tt)
      call smd_data_create_get(namespace,tag,na,MT_STR,i_tt)
      tag = "param:charge"
      call smd_data_create_get(namespace,tag,na,MT_DBL,i_q)
      tag = "param:mass"
      call smd_data_create_get(namespace,tag,na,MT_DBL,i_m)
      tag = "param:resname"
c      call smd_data_create_get(namespace,tag,16*na,MT_BYTE,i_tr)
      call smd_data_create_get(namespace,tag,na,MT_STR,i_tr)

c
c     fill in names and coordinates from parameter file
c     -------------------------------------------------

      call smd_paramfile_read_atom(na,
     +                         byte_mb(i_ta),
     +                         byte_mb(i_tt),
     +                         byte_mb(i_tr),
     +                         dbl_mb(i_m),
     +                         dbl_mb(i_q))

      return
      end

      subroutine smd_param_ntot(namespace,ntot)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) namespace
      integer ntot
c
      character*72 tag
      character*30 pname
      logical result

      pname = "smd_param_ntot"
      tag = "param:charge"
      call smd_data_get_dim(namespace,tag,ntot,result)
      if(.not. result)  
     >  call errquit(
     >       pname//'error getting size for '//tag,0, 0)
      ntot = ntot

      return
      end

c $Id$
