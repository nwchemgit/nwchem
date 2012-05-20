      subroutine smd_cutoff_init()
      implicit none
#include "rtdb.fh"
#include "srtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"

c
      character*32 tag
      character*32 pname
      integer srtdb,rtdb
      double precision r

      pname = "smd_cutoff_init"

      call smd_rtdb_get_handle(rtdb)
      call smd_rtdb_get_handle(srtdb)
      tag="smd:rcut"
      if (.not.rtdb_get(rtdb,tag,mt_dbl,1,r)) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)
      if (.not.srtdb_put(srtdb,"rcut",mt_dbl,1,r)) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)

      tag="smd:verlet" 
      if (.not.rtdb_get(rtdb,tag,mt_dbl,1,r)) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)
      if (.not.srtdb_put(rtdb,"rcut_verlet",mt_dbl,1,r)) 
     >      call errquit(pname//'failed to get'//tag,0,
     >       RTDB_ERR)


      end 
c $Id$
