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
      write(*,*) "number of atoms in parameter space",na
      if(na.le.0)
     >  call errquit(
     >       pname//'no parameters fond',0, 0)

c
      call smd_namespace_create(namespace)
c
c     create parameter data structures
c     ---------------------------
      call smd_data_create(namespace,"param:atomname",16*na,MT_BYTE)
      call smd_data_create(namespace,"param:typename",16*na,MT_BYTE)
      call smd_data_create(namespace,"param:charge",na,MT_DBL)
      call smd_data_create(namespace,"param:mass",na,MT_DBL)
      call smd_data_create(namespace,"param:resname",16*na,MT_BYTE)

c
c     fill in names and coordinates from parameter file
c     -------------------------------------------------
      tag = "param:atomname"
      call smd_data_get_index(namespace,tag,i_ta,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "param:typename"
      call smd_data_get_index(namespace,tag,i_tt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      tag = "param:resname"
      call smd_data_get_index(namespace,tag,i_tr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      tag = "param:mass"
      call smd_data_get_index(namespace,tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

      tag = "param:charge"
      call smd_data_get_index(namespace,tag,i_q,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)


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
      call smd_data_get_size(namespace,tag,ntot,result)
      if(.not. result)  
     >  call errquit(
     >       pname//'error getting size for '//tag,0, 0)
      ntot = ntot

      return
      end

c $Id$
