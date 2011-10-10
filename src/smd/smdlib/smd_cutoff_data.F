      subroutine smd_cutoff_get_rcut(arcut)
      implicit none
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      double precision arcut
c
      integer i
      character*(smd_tag_size) tag
      logical result
c
c      arcut = rcut
      tag = "rcut"
      call smd_get_ind(tag,i,result)
      if(.not. result)
     >  call errquit(
     >       'error getting index for'//tag,0, 0)
      arcut = dbl_mb(i)

      end 

      subroutine smd_cutoff_get_rcut_verlet(arcut)
      implicit none
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      double precision arcut
c
      integer i
      character*(smd_tag_size) tag
      logical result
c
      tag = "rcut_verlet"
      call smd_get_ind(tag,i,result)
      if(.not. result)
     >  call errquit(
     >       'error getting index for'//tag,0, 0)
      arcut = dbl_mb(i)

      end 

      subroutine smd_cutoff_init()
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"

c
      character*(smd_tag_size) namespace,tag
      character*32 pname
      logical result
      integer i
c
      double precision rcut,rcut_verlet
      pname = "smd_cutoff_init"

      tag="smd:rcut"
      call smd_rtdb_get_dbl(tag,1,rcut,result)
      if (.not.result) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)

      tag="smd:verlet" 
      call smd_rtdb_get_dbl(tag,1,rcut_verlet,result)
      if (.not.result) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)

      tag = "cutoff"
      call smd_system_get_component(namespace,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_namespace_create(namespace)
      tag = "rcut"
      call smd_data_create_get(namespace,tag,1,MT_DBL,i)
      dbl_mb(i) = rcut
      tag = "rcut_verlet"
      call smd_data_create_get(namespace,tag,1,MT_DBL,i)
      dbl_mb(i) = rcut_verlet


      end 
c $Id$
