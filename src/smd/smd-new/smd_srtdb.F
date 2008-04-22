      subroutine smd_srtdb_init()
      implicit none
#include "rtdb.fh"
#include "srtdb.fh"
#include "smd_srtdb_data.fh"
#include "mafdecls.fh"
c

      if (.not. srtdb_open('smd.main', 'empty', smd_srtdb)) 
     *     call errquit('srtdb_open',0,0)

      end 


      subroutine smd_srtdb_get_handle(srtdb)
      implicit none
#include "rtdb.fh"
#include "smd_srtdb_data.fh"
#include "mafdecls.fh"
c
      integer srtdb
c
      srtdb = smd_srtdb

      end 
