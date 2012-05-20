      subroutine smd_rtdb_get_handle(rtdb)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      integer rtdb
c
      rtdb= smd_rtdb

      end 

      subroutine smd_rtdb_get_istart(istart)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      integer istart
c
      istart= smd_istart

      end 

      subroutine smd_rtdb_get_iend(iend)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      integer iend
c
      iend= smd_iend

      end 

      subroutine smd_rtdb_get_nproc(nproc)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      integer nproc
c
      nproc= smd_nproc

      end 

      subroutine smd_rtdb_get_name(rtdb_name)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      character*(*) rtdb_name
c
      rtdb_name= smd_rtdb_name

      end 

      subroutine smd_rtdb_set_handle(rtdb)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      integer rtdb
c
      smd_rtdb = rtdb

      end 

      subroutine smd_rtdb_get_operiodic(operiodic)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      logical operiodic
c
      double precision latt(3,3)
      character*32 pname
      character*80 tag
      double precision a(3)
      integer i
c
      pname = "smd_lat_rtdb_read"
c
c      write(*,*) "in "//pname
c
      tag="smd:operiodic"
      if (rtdb_get(smd_rtdb,tag,mt_log,1,operiodic)) 
     >    return

      operiodic = .true.
      tag="smd:lat_a"
      if (.not.rtdb_get(smd_rtdb,tag,mt_dbl,3,a(1))) 
     >      operiodic=.false.


      end 

      subroutine smd_rtdb_get_paramfile(filename,result)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      character*(*) filename
      logical result
c
      
      result = .true.
      if(.not.rtdb_cget(smd_rtdb,'smd:paramfile',1,filename))
     >   result = .false.
      end 

      subroutine smd_rtdb_get_coordfile(filename,result)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "global.fh"

c
      character*(*) filename
      logical result
c
      character*30 pname

      pname = "smd_rtdb_get_coordfile"
      
      result = .true.
      if(.not.rtdb_cget(smd_rtdb,'smd:coordfile',1,filename))
     >   result = .false.

      end 

      subroutine smd_rtdb_init(rtdb)
      implicit none
#include "rtdb.fh"
#include "srtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"

      integer rtdb
c
      logical ignore
      smd_rtdb = rtdb
      if(.not. rtdb_getfname(rtdb, smd_rtdb_name)) call
     *        errquit('rtdb_getfname failed',0,0)

      smd_istart=0
      smd_iend  =0
      smd_nproc =0
      ignore = rtdb_get(rtdb,"smd:istart",mt_int,1,smd_istart)
      ignore = rtdb_get(rtdb,"smd:iend",mt_int,1,smd_iend)
      ignore = rtdb_get(rtdb,"smd:nproc",mt_int,1,smd_nproc)

      if (.not. srtdb_open('smd.main', 'empty', smd_srtdb)) 
     *     call errquit('srtdb_open',0,0)

      if (.not. srtdb_put(rtdb, 'smd_istart', MT_INT, 1, smd_istart))
     $     call errquit('srtdb_put',0,0)

      if (.not. srtdb_get(rtdb, 'smd_istart', MT_INT, 1, smd_istart))
     $     call errquit('srtdb_put',0,0)

      write(*,*) "found istart in srtdb",smd_istart
      end 
c $Id$
