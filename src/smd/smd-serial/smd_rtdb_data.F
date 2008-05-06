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

      subroutine smd_rtdb_get_veloc_input(filename,result)
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

      pname = "smd_rtdb_get_veloc:input"
      
      result = .true.
      if(.not.rtdb_cget(smd_rtdb,'smd:veloc:input',1,filename))
     >   result = .false.

      end 

      subroutine smd_rtdb_init(rtdb)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "smd_const_data.fh"

      integer rtdb
c
      character*(smd_string_size) namespace
      character*(smd_string_size) tag
      logical ignore,result
      character*30 pname

      pname = "smd_rtdb_init"
      smd_rtdb = rtdb
      if(.not. rtdb_getfname(rtdb, smd_rtdb_name)) call
     *        errquit('rtdb_getfname failed',0,0)

      smd_istart=0
      smd_iend  =0
      smd_nproc =0
      ignore = rtdb_get(rtdb,"smd:istart",mt_int,1,smd_istart)
      ignore = rtdb_get(rtdb,"smd:iend",mt_int,1,smd_iend)
      ignore = rtdb_get(rtdb,"smd:nproc",mt_int,1,smd_nproc)

      tag = "rtdb"
      call smd_system_get_component(namespace,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_namespace_create(namespace)

      call smd_data_namespace_rtdb_get(rtdb,namespace,"smd")
      call smd_data_namespace_print(namespace,6)
      end 
