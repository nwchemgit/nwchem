            block data initial_smd_cutoff_data
      implicit none
#include "smd_cutoff_data.fh"
c     
      data ostatus /.false./
      data rcut /0.0d0/
      data rcut_verlet /0.0d0/


      end

      subroutine smd_cutoff_get_rcut(arcut)
      implicit none
#include "smd_cutoff_data.fh"

      double precision arcut

      arcut = rcut

      end 

      subroutine smd_cutoff_set_rcut(arcut)
      implicit none
#include "smd_cutoff_data.fh"

      double precision arcut

      rcut = arcut

      end 

      subroutine smd_cutoff_get_rcut_verlet(arcut_verlet)
      implicit none
#include "smd_cutoff_data.fh"

      double precision arcut_verlet

      arcut_verlet = rcut_verlet

      end 

      subroutine smd_cutoff_set_rcut_verlet(arcut_verlet)
      implicit none
#include "smd_cutoff_data.fh"

      double precision arcut_verlet

      rcut_verlet = arcut_verlet

      end 

      subroutine smd_cutoff_print(un)
      implicit none
#include "smd_cutoff_data.fh"

      integer un

      write(un,*) "cutoff data"
      write(un,*) "cutoff radius",rcut
      write(un,*) "verlet cutoff radius",rcut_verlet

      end 

      subroutine smd_cutoff_init()
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_cutoff_data.fh"

c
      character*32 tag
      character*32 pname
      integer rtdb

      pname = "smd_cutoff_init"

      call smd_rtdb_get_handle(rtdb)
      tag="smd:rcut"
        if (.not.rtdb_get(rtdb,tag,mt_dbl,1,rcut)) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)

      tag="smd:verlet" 
        if (.not.rtdb_get(rtdb,tag,mt_dbl,1,rcut_verlet)) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)


      end 
c $Id$
