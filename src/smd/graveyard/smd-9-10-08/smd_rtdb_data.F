      subroutine smd_rtdb_close(action)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      character*(*) action
c
      integer rtdb
      logical result
      character*30 pname

      pname = "smd_rtdb_close"
      call smd_rtdb_get_handle(rtdb)
      if(rtdb.lt.0) call errquit(pname//"no rtdb to close")
      if(.not.rtdb_close(rtdb,action)) 
     >     call errquit(pname//" while rtdb_close",0,0)
      call smd_rtdb_set_handle(-1000)

      end 

      subroutine smd_rtdb_open(fname,action)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      character*(*) fname
      character*(*) action
c
      integer rtdb
      logical result
      character*30 pname

      pname = "smd_rtdb_close"
      call smd_rtdb_get_handle(rtdb)
      if(rtdb.ge.0) call errquit(pname//"close rtdb first")
      if(.not.rtdb_open(fname,action,rtdb)) 
     >     call errquit(pname//" while rtdb_open",0,0)
      call smd_rtdb_set_handle(rtdb)

      end 

      subroutine smd_rtdb_get_handle(rtdb)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      integer rtdb
c
      logical result
      character*30 pname

      pname = "smd_rtdb_get_handle"
      call smd_rtdb_get_int("smd:rtdb_handle",1,rtdb,result)
      if(.not.result) call errquit(pname//"no rtdb component",0,0)

      end 

      subroutine smd_rtdb_set_handle(rtdb)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      integer rtdb
c
      logical result
      character*30 pname

      pname = "smd_rtdb_set_handle"
      call smd_rtdb_put_int("smd:rtdb_handle",1,rtdb,result)
      if(.not.result) call errquit(pname//"no rtdb component",0,0)

      end 

      subroutine smd_rtdb_get_istart(istart)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      integer istart
c
      logical result
      character*30 pname
c
      pname = "smd_rtdb_get_istart"
      call smd_rtdb_get_int("smd:fragment_istart",1,istart,result)
      if(.not.result) call errquit(pname//"no rtdb component",0,0)

      end 

      subroutine smd_rtdb_get_iend(iend)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      integer iend
c
      logical result
      character*30 pname
c
      pname = "smd_rtdb_get_iend"
      call smd_rtdb_get_int("smd:fragment_iend",1,iend,result)
      if(.not.result) call errquit(pname//"no rtdb component",0,0)

      end 

      subroutine smd_rtdb_get_nproc(nproc)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
c
      integer nproc
c
      logical result
      character*30 pname
c
      pname = "smd_rtdb_get_nproc"
      call smd_rtdb_get_int("smd:fragment_nproc",1,nproc,result)
      if(.not.result) call errquit(pname//"no rtdb component",0,0)

      end 

      subroutine smd_rtdb_get_name(rtdb_name)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      character*(*) rtdb_name
c
      integer rtdb
      character*30 pname
      logical result

      pname = "smd_rtdb_get_name"
c      call smd_rtdb_get_handle(rtdb)
c      if(.not. rtdb_getfname(rtdb, rtdb_name)) call
c     *        errquit(pname//'rtdb_getfname failed',0,0)
      call smd_rtdb_get_string("smd:rtdb_name",1,rtdb_name,result)
      if(.not.result) call errquit(pname//"no rtdb component",0,0)


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

      subroutine smd_rtdb_init(parallel,rtdb)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "smd_const_data.fh"

      logical parallel
      integer rtdb
c
      character*(smd_string_size) namespace
      character*(smd_string_size) tag
      logical ignore,result,oldmode
      character*30 pname

      pname = "smd_rtdb_init"
      oldmode = rtdb_parallel(parallel)
      smd_rtdb = rtdb
c      if(.not. rtdb_getfname(rtdb, smd_rtdb_name)) call
c     *        errquit('rtdb_getfname failed',0,0)

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
      oldmode = rtdb_parallel(oldmode)
      end 

      subroutine smd_rtdb_get_log(dname,ndim,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
      character*(smd_string_size) sname
      character*(*) dname
      integer ndim
      logical datum(ndim)
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_rtdb_get_log"

      call smd_system_get_component(sname,"rtdb",result)
      if(.not.result) call errquit(pname//"no rtdb component",0,0)
      call smd_data_get_log(sname,dname,ndim,datum,result)
      return
      end

      subroutine smd_rtdb_get_int(dname,ndim,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
      character*(smd_string_size) sname
      character*(*) dname
      integer ndim
      integer datum(ndim)
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_rtdb_get_int"

      call smd_system_get_component(sname,"rtdb",result)
      if(.not.result) call errquit(pname//"no rtdb component",0,0)
      call smd_data_get_int(sname,dname,ndim,datum,result)
      return
      end

      subroutine smd_rtdb_put_int(dname,ndim,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
      character*(smd_string_size) sname
      character*(*) dname
      integer ndim
      integer datum(ndim)
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_rtdb_put_int"

      call smd_system_get_component(sname,"rtdb",result)
      if(.not.result) return
      call smd_data_put_int(sname,dname,ndim,datum)
      return
      end

      subroutine smd_rtdb_get_dbl(dname,ndim,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
      character*(smd_string_size) sname
      character*(*) dname
      integer ndim
      double precision datum(ndim)
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_rtdb_get_dbl"
     
      call smd_system_get_component(sname,"rtdb",result)
      if(.not.result) call errquit(pname//"no rtdb component",0,0)
      call smd_data_get_dbl(sname,dname,ndim,datum,result)
      return
      end

      subroutine smd_rtdb_get_string(dname,ndim,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
#include "inp.fh"
      character*(*) dname
      integer ndim
      integer dtype
      character*(*) datum(ndim)
c
      integer ind
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1
      integer chunk
      logical result
      character*(smd_string_size) sname

      pname = "smd_rtdb_get_string"

      call smd_system_get_component(sname,"rtdb",result)
      if(.not.result) call errquit(pname//"no rtdb component",0,0)
      call smd_data_get_string(sname,dname,ndim,datum,result)
      return
      end

c $Id$
