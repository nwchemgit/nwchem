      subroutine smd_rtdb_init(parallel,rtdb)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "smd_const_data.fh"

      logical parallel
      integer rtdb
c
      character*(smd_string_size) namespace
      character*(smd_string_size) tag
      logical ignore,result,oldmode
      character*256 big_buffer
      character*(smd_string_size) buffer1
      logical oscratch(smd_rtdb_nmax)
      character*30 pname
      integer i
      integer istart,iend

      pname = "smd_rtdb_init"
      oldmode = rtdb_parallel(parallel)
      smd_rtdb = rtdb
c
      smd_istart=0
      smd_iend  =0
      smd_nproc =0
      ignore = rtdb_get(rtdb,"smd:istart",mt_int,1,smd_istart)
      ignore = rtdb_get(rtdb,"smd:iend",mt_int,1,smd_iend)
      ignore = rtdb_get(rtdb,"smd:nproc",mt_int,1,smd_nproc)

c
c     create rtdb namespace
c     ---------------------
      tag = "rtdb"
      call smd_system_get_component(namespace,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
      call smd_namespace_create(namespace)
c
c     fetch and store rtdb name
c     -------------------------
      if(.not. rtdb_getfname(rtdb, big_buffer)) call
     *        errquit('rtdb_getfname failed',0,0)
      istart = 0
   10 if (inp_strtok(big_buffer, '/', istart, iend)) then
      buffer1 = big_buffer(istart:iend)
      goto 10
      endif
c
      tag  = "smd:rtdb_name"
      call smd_data_create(namespace,tag,smd_rtdb_nmax,
     >                         MT_STR)
      call smd_data_put_string1(namespace,tag,1,1,
     >                              buffer1)
c
c     store prefix
c     ------------
      call util_file_prefix_get(buffer1)
      tag  = "smd:prefix"
      call smd_data_create(namespace,tag,smd_rtdb_nmax,
     >                         MT_STR)
      call smd_data_put_string1(namespace,tag,1,1,
     >                         buffer1)
c
c     store location of rtdb (perm by default)
c     ---------------------------------------
      tag  = "smd:rtdb_scratch"
      call smd_data_create(namespace,tag,smd_rtdb_nmax,
     >                         MT_LOG)
      do i=1,smd_rtdb_nmax
        oscratch(i)=.false.
      end do
      call smd_data_put_log1(namespace,tag,1,1,
     >                         .false.)
c
c     store rtdb pointer
c     ------------------
      tag  = "smd:rtdb_pointer"
      call smd_data_create(namespace,tag,1,
     >                         MT_INT)
      call smd_data_put_int1(namespace,tag,1,1,
     >                         1)
c
c     store rtdb handle
c     ------------------
      tag  = "smd:rtdb_handle"
      call smd_data_create(namespace,tag,1,
     >                         MT_INT)
      call smd_data_put_int1(namespace,tag,1,1,
     >                         rtdb)
c
c     digest some of rtdb contents
c     ----------------------------
      call smd_data_namespace_rtdb_get(rtdb,namespace,"smd")
      call smd_data_namespace_rtdb_get(rtdb,namespace,"geom:source")
      call smd_data_namespace_rtdb_get(rtdb,namespace,"task:theory")
      call smd_data_namespace_rtdb_get(rtdb,namespace,"task:operation")
c
c     restore rtdb access mode
c     -----------------------
      oldmode = rtdb_parallel(oldmode)
c    
      return
      end 

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

      subroutine smd_rtdb_open(rtdb)
c     opens rtdb whose name was already stored
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
c
      integer rtdb
c
      logical result
      character*30 pname
      character*(smd_path_size) rtdb_name
      character*(smd_string_size) prefix
      logical oscratch

      pname = "smd_rtdb_open"
c     check if other rtdb is open
c     ---------------------------
      call smd_rtdb_get_handle(rtdb)
      if(rtdb.ge.0) call errquit(pname//"close rtdb first",0,0)
c     get new rtdb name
c     -----------------
      call smd_rtdb_get_name_full(rtdb_name)
c      note that mode is to old so the actual file should be there!
      if(.not.rtdb_open(rtdb_name,"old",rtdb)) 
     >     call errquit(pname//" cannot find rtdb "//rtdb_name,0,0)
      call smd_rtdb_set_handle(rtdb)
      call smd_rtdb_get_prefix(prefix)
      call util_file_prefix_set(prefix)

      end 

      subroutine smd_rtdb_push_open_clone(taskid,oscratch,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "rtdb.fh"
      integer taskid
      logical oscratch
      integer rtdb
c     
      character*30 pname
c
      pname = "smd_rtdb_push_open_clone"
      call smd_rtdb_push_clone(taskid,oscratch)
      call smd_rtdb_open(rtdb)
      call ga_sync()

      return
      end

      subroutine smd_rtdb_push_clone(taskid,oscratch)
c     rtdb cloning
c     1. make a numbered name (taskid)
c     2. copy file
c
c     taskid[IN] task number
c     rtdb[OUT]  rtdb handle
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "rtdb.fh"
      integer taskid
      logical oscratch
c     
      character*(smd_string_size)  sname
      character*(smd_string_size)  rtdb_name0,rtdb_name
      character*(smd_path_size)  rtdb_name_full0
      character*(smd_path_size)  rtdb_name_full
      integer i1,i2
      logical result
      logical master
      character*(smd_string_size) prefix, prefix0
      integer rtdb
      integer n,nmax
      character*30 pname
c
      pname = "smd_rtdb_push_clone1"
c      
      call ga_sync()
      master = ga_nodeid().eq.0
c     get full path to current rtdb 
      call smd_rtdb_get_name_full(rtdb_name_full0)
c
c     create new prefix name
c     ----------------------
      call smd_rtdb_get_prefix(prefix0)
      i1 = inp_strlen(prefix0)
      write(prefix,'(A,I4.4)')
     >    prefix0(1:i1),taskid
      i1 = inp_strlen(prefix)
      write(rtdb_name,'(A,A3)')
     >    prefix(1:i1),".db"
c
c     push names in 
c     -------------
      call util_flush(6)
      call smd_rtdb_push_name(rtdb_name,prefix,oscratch)
c
c     get full path to new rtdb 
c     ----------------------------
      call smd_rtdb_get_name_full(rtdb_name_full)
c
c     actually copy the file
c     ----------------------
      if(master) then
        inquire(file=rtdb_name_full0,exist=result)
        if(.not.result) then
          call errquit("cannot find rtdb"//
     >                 rtdb_name_full0,
     >                 0,0)
        end if
        call util_file_copy(rtdb_name_full0,rtdb_name_full)
      end if
c
      call ga_sync()

      return
      end

      subroutine smd_rtdb_pop()
c
      implicit none
#include "errquit.fh"
      character*30 pname
c
      pname = "smd_rtdb_pop"
c      
      call smd_rtdb_close("delete")
      call smd_rtdb_pop_name()
      return
      end

      subroutine smd_rtdb_get_handle(rtdb)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
c
      integer rtdb
c
      logical result
      character*30 pname
      character*(smd_string_size) tag,sname

      pname = "smd_rtdb_get_handle"
      tag = "rtdb"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "smd:rtdb_handle"
      call smd_data_get_int1(sname,tag,1,1,rtdb)

      end 

      subroutine smd_rtdb_set_handle(rtdb)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
c
      integer rtdb
c
      logical result
      character*30 pname
      character*(smd_string_size) tag,sname

      pname = "smd_rtdb_get_handle"
      tag = "rtdb"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "smd:rtdb_handle"
      call smd_data_put_int1(sname,tag,1,1,rtdb)

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

      subroutine smd_rtdb_get_prefix(aname)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      character*(*) aname
c
      character*(smd_string_size)  buffer
      character*(smd_string_size)  tag,sname
      integer rtdb
      character*30 pname
      integer n
      logical result

      pname = "smd_rtdb_get_prefix"
      tag = "rtdb"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_data_get_int1(sname,"smd:rtdb_pointer",1,1,
     >                         n)

      call smd_data_get_string1(sname,"smd:prefix",n,n,
     >                          buffer)


      aname = buffer

      end 

      subroutine smd_rtdb_get_name(rtdb_name)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      character*(*) rtdb_name
c
      character*(smd_string_size)  buffer
      character*(smd_string_size)  tag,sname
      integer rtdb
      character*30 pname
      integer n
      logical result

      pname = "smd_rtdb_get_name"
      tag = "rtdb"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_data_get_int1(sname,"smd:rtdb_pointer",1,1,
     >                         n)

      call smd_data_get_string1(sname,"smd:rtdb_name",n,n,
     >                          buffer)


      rtdb_name = buffer

      end 

      subroutine smd_rtdb_get_name_full(rtdb_name)
c     obtain current rtdb name with fully resolved path
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      character*(*) rtdb_name
c
      character*(smd_string_size)  tag,sname
      character*30 pname
      integer n
      logical result
      logical oscratch

      pname = "smd_rtdb_get_name_full"
      tag = "rtdb"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "smd:rtdb_pointer"
      call smd_data_get_int1(sname,tag,1,1,n)

      tag = "smd:rtdb_name"
      call smd_data_get_string1(sname,tag,n,n,rtdb_name)

      tag  = "smd:rtdb_scratch"
      call smd_data_get_log1(sname,tag,n,n,oscratch)

      call util_file_name_resolve(rtdb_name,oscratch)

      end 

      subroutine smd_rtdb_push_name(rtdb_name,prefix,oscratch)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      character*(*) rtdb_name,prefix
      logical oscratch
c
      character*(smd_string_size)  buffer(smd_rtdb_nmax)
      character*(smd_string_size)  tag,sname
      integer rtdb
      character*30 pname
      integer n
      logical result

      pname = "smd_rtdb_push_name"
      tag = "rtdb"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      tag = "smd:rtdb_pointer"
      call smd_data_get_int1(sname,tag,1,1, n)
      if(.not.result) call errquit(pname//"no component "//tag,0,0)
c
      n=n+1
      if(n.gt.smd_rtdb_nmax) 
     > call errquit(pname//"out of stack",0,0)
      tag = "smd:rtdb_pointer"
      call smd_data_put_int1(sname,tag,1,1,n)

      tag = "smd:rtdb_name"
      call smd_data_put_string1(sname,tag,n,n,rtdb_name)

      tag = "smd:prefix"
      call smd_data_put_string1(sname,tag,n,n,prefix)

      tag  = "smd:rtdb_scratch"
      call smd_data_put_log1(sname,tag,n,n,oscratch)

      return
      end 

      subroutine smd_rtdb_pop_name()
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      character*(smd_string_size)  tag,sname
      integer rtdb
      character*30 pname
      integer n
      logical result

      pname = "smd_rtdb_pop_name"
      tag = "rtdb"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "smd:rtdb_pointer"
      call smd_data_get_int1(sname,tag,1,1, n)
      if(.not.result) call errquit(pname//"no component "//tag,0,0)
      if(n-1.le.0) 
     > call errquit(pname//"out of stack",0,0)
c
      tag = "smd:rtdb_name"
      call smd_data_put_string1(sname,tag,n,n," ")

      tag = "smd:prefix"
      call smd_data_put_string1(sname,tag,n,n," ")
c
      n=n-1
      tag = "smd:rtdb_pointer"
      call smd_data_put_int1(sname,tag,1,1,n)

      return
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
