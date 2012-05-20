      subroutine smd_temper_get_ke_target(k)
      implicit none
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      double precision k
c
      integer degfree, natms,ntcons,i
      character*(smd_tag_size) tag
      double precision t
      logical result
c
      tag = "target_temperature"
      call smd_get_ind(tag,i,result)
      if(.not. result)
     >  call errquit(
     >       'error getting index for'//tag,0, 0)
      t = dbl_mb(i)


      call smd_atom_ntot(natms)
      call smd_system_ntcons(ntcons)

      degfree = dble(3*(natms)-3-ntcons)
      k       = degfree*t*boltzmann*0.5


      end 

      subroutine smd_temper_init()
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
c
      character*(smd_tag_size) namespace,tag
      character*32 pname
      integer natms,ntcons
      logical result
      integer i
      double precision t
c
      pname = "smd_temper_init"

      tag="smd:temp_target"
      call smd_rtdb_get_dbl(tag,1,t,result)
      if (.not.result) then
        call smd_system_unset_component("temperature")
         return
      end if

      tag = "temperature"
      call smd_system_get_component(namespace,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      call smd_namespace_create(namespace)
      tag = "target_temperature"
      call smd_data_create_get(namespace,tag,1,MT_DBL,i)
      dbl_mb(i) = t

      tag = "current_temperature"
      call smd_data_create_get(namespace,tag,1,MT_DBL,i)
      dbl_mb(i) = -1.0d0

      return

      end 

      subroutine smd_temper_ke_compute(t,ke)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
c
      double precision t,ke
      integer degfree
      integer natms,ntcons
c
      call smd_atom_ntot(natms)
      call smd_system_ntcons(ntcons)

      degfree = dble(3*(natms)-3-ntcons)
      ke      = degfree*t*boltzmann*0.5


      end 
c $Id$
