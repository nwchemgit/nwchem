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

      subroutine smd_srtdb_get_paramfile(filename,result)
      implicit none
#include "srtdb.fh"
#include "smd_srtdb_data.fh"
#include "mafdecls.fh"
c
      character*(*) filename
      logical result
c
      
      result = .true.
      if(.not.srtdb_cget(smd_srtdb,'smd:paramfile',1,filename))
     >   result = .false.
      end 

      subroutine smd_srtdb_get_coordfile(filename,result)
      implicit none
#include "srtdb.fh"
#include "smd_srtdb_data.fh"
#include "mafdecls.fh"
#include "global.fh"

c
      character*(*) filename
      logical result
c
      character*30 pname

      pname = "smd_srtdb_get_coordfile"
      
      result = .true.
      if(.not.srtdb_cget(smd_srtdb,'smd:coordfile',1,filename))
     >   result = .false.

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
c $Id$
