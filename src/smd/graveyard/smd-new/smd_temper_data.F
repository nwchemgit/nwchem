            block data initial_smd_temper_data
      implicit none
#include "smd_temper_data.fh"
c     
      data otemp /.false./
      data temp_target /0.0d0/

      end

      subroutine smd_temper_get_temp_target(t)
      implicit none
#include "smd_temper_data.fh"

      double precision t

      t = temp_target

      end 

      subroutine smd_temper_get_ke_target(k)
      implicit none
#include "smd_temper_data.fh"

      double precision k

      k = ke_target

      end 

      subroutine smd_temper_init()
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_temper_data.fh"
#include "smd_system.fh"
#include "smd_const_data.fh"
c
      character*32 tag
      character*32 pname
      integer rtdb
      integer natms,ntcons

      pname = "smd_temper_init"

      call smd_rtdb_get_handle(rtdb)
      tag="smd:temp_target"
        if (.not.rtdb_get(rtdb,tag,mt_dbl,1,temp_target)) 
     >      temp_target = 0.0d0

      if(.not.smd_system_ntot(natms))
     >      call errquit(pname//'failed to get ntot',0,
     >       RTDB_ERR)

      if(.not.smd_system_ntcons(ntcons))
     >      call errquit(pname//'failed to get ntcons',0,
     >       RTDB_ERR)

      degfree = dble(3*(natms)-3-ntcons)
      ke_target = degfree*temp_target*boltzmann*0.5

      end 

      subroutine smd_temper_ke_compute(t,ke)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_temper_data.fh"
#include "smd_system.fh"
#include "smd_const_data.fh"
c
      double precision t,ke

      ke = degfree*t*boltzmann*0.5

      end 
c $Id$
