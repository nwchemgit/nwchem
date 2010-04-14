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
      character*72 buffer(smd_rtdb_nmax)
      logical oscratch(smd_rtdb_nmax)
      character*30 pname
      integer i
      integer istart,iend

      pname = "smd_rtdb_init"
      write(*,*) "in "//pname
      oldmode = rtdb_parallel(parallel)
      smd_rtdb = rtdb
c
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


c
c     fetch and store rtdb name
c     -------------------------
      if(.not. rtdb_getfname(rtdb, big_buffer)) call
     *        errquit('rtdb_getfname failed',0,0)
      do i=1,smd_rtdb_nmax
        buffer(i) = " "
        oscratch(i)=.false.
      end do
      istart = 0
   10 if (inp_strtok(big_buffer, '/', istart, iend)) then
      buffer(1) = big_buffer(istart:iend)
      goto 10
      endif
c
c      call smd_data_put_string(namespace,"smd:rtdb_name",smd_rtdb_nmax,
c     >                         buffer)
      tag  = "smd:rtdb_name"
      call smd_data_create(namespace,tag,smd_rtdb_nmax,
     >                         MT_STR)
      call smd_data_put_string1(namespace,tag,1,1,
     >                              buffer(1))
c     -------------------------
      call util_file_prefix_get(buffer(1))
      tag  = "smd:prefix"
      call smd_data_create(namespace,tag,smd_rtdb_nmax,
     >                         MT_STR)
      call smd_data_put_string1(namespace,tag,1,1,
     >                         buffer(1))

      call smd_data_put_int(namespace,"smd:rtdb_pointer",1,
     >                         1)

      call smd_data_namespace_print(namespace,50)
      call smd_data_namespace_rtdb_get(rtdb,namespace,"smd")
      call smd_data_namespace_rtdb_get(rtdb,namespace,"geom:source")
      oldmode = rtdb_parallel(oldmode)

      stop
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

      pname = "smd_rtdb_open"
      call smd_rtdb_get_handle(rtdb)
      if(rtdb.ge.0) call errquit(pname//"close rtdb first")
      if(.not.rtdb_open(fname,action,rtdb)) 
     >     call errquit(pname//" while rtdb_open",0,0)
      write(*,*) "rtdb handle is",rtdb
      call smd_rtdb_set_handle(rtdb)

      end 

      subroutine smd_rtdb_push_clone(taskid)
c     rtdb cloning
c     1. make a numbered name (taskid)
c     2. copy file
c     3. open rtdb
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
c     
      character*(smd_string_size)  sname
      character*(smd_string_size)  rtdb_name0,rtdb_name
      character*256  rtdb_name_full0
      character*256  rtdb_name_full
      integer i1,i2
      logical result
      logical master
      character*(smd_string_size) prefix, prefix0
      integer rtdb
      integer n,nmax
      character*30 pname
c
      pname = "smd_rtdb_push_clone"
c      
      write(*,*) "in "//pname
      call ga_sync()
      master = ga_nodeid().eq.0
c
c     create new prefix name
c     ----------------------
      call smd_rtdb_get_prefix1(prefix0)
      i1 = inp_strlen(prefix0)
      write(prefix,'(A,I4.4)')
     >    prefix0(1:i1),taskid
      call util_file_prefix_set(prefix)
c
c     create new rtdb name
c     --------------------
      call util_file_prefix("db",rtdb_name)
      rtdb_name_full=rtdb_name
      write(*,*) "new db name ",rtdb_name_full
      call util_file_name_resolve(rtdb_name_full,.false.)
      write(*,*) "new db name ",rtdb_name_full
c
c     get absolute name of the original database
c     ------------------------------------------
      call smd_rtdb_get_name1(rtdb_name_full0)
      write(*,*) "old db name ",rtdb_name_full0
      call util_file_name_resolve(rtdb_name_full0,.false.)
      write(*,*) "old db name ",rtdb_name_full0
c
      if(master) then
        inquire(file=rtdb_name_full0,exist=result)
        if(.not.result) then
          call errquit("cannot find rtdb"//
     >                 rtdb_name_full0,
     >                 0,0)
        end if
        call util_file_copy(rtdb_name_full0,rtdb_name_full)
      end if

      call smd_rtdb_close("keep")
      call smd_rtdb_open(rtdb_name_full,"old")
      call smd_rtdb_push_name(prefix,rtdb_name)

      call ga_sync()
      write(*,*) "out "//pname

      return
      end

      subroutine smd_rtdb_open_clone(taskid,rtdb)
c     rtdb cloning
c     1. make a numbered name (taskid)
c     2. copy file
c     3. open rtdb
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
      integer taskid,rtdb
c     
      character*(smd_string_size)  buffer
      character*256  rtdb_name
      character*256  rtdb_name_group
      integer i1,i2
      logical result
      logical master
      character*255 prefix, prefix0
      character*30 pname
c
      pname = "smd_task_rtdb_clone"
c
      call ga_sync()
      master = ga_nodeid().eq.0
c
      call smd_rtdb_get_prefix(prefix0)
      i1 = inp_strlen(prefix0)
      write(prefix,'(A,I4.4)')
     >    prefix0(1:i1),taskid
      call util_file_prefix_set(prefix)
      call util_message("prefix is set to "//prefix)
c     get absolute name of the original database
      call smd_rtdb_get_name(rtdb_name)
      write(*,*) "rtdb_name",rtdb_name
      call util_file_name_resolve(rtdb_name,.false.)
      write(*,*) "resolved rtdb_name",rtdb_name
c     create name for group database in scratch
c      write(buffer,'(I4.4,A3)') taskid,".db"
      write(buffer,'(A2)') "db"
      call util_file_name(buffer,.true.,.false.,rtdb_name_group)
      write(*,*) "resolved group_rtdb_name",rtdb_name_group
      call util_flush(6)
c
      if(master) then
        inquire(file=rtdb_name,exist=result)
        if(.not.result) then
          call errquit("cannot find rtdb"//
     >                 rtdb_name,
     >                 0,0)
        end if
        call util_file_copy(rtdb_name,rtdb_name_group)
      end if

      call smd_rtdb_open(rtdb_name_group,"old")
      call smd_rtdb_get_handle(rtdb)
      write(*,*) "got new rtdb handle",rtdb
      call ga_sync()

      return
      end

      subroutine smd_rtdb_open_master(rtdb)
c     rtdb cloning
c     1. make a numbered name (taskid)
c     2. copy file
c     3. open rtdb
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
      integer taskid,rtdb
c     
      character*(smd_string_size)  buffer
      character*256  rtdb_name
      character*256  rtdb_name_group
      integer i1,i2
      logical result
      logical master
      character*255 prefix, prefix0
      character*30 pname
c
      pname = "smd_task_rtdb_clone"
c
      call ga_sync()
      master = ga_nodeid().eq.0
c
c     get absolute name of the original database
      call smd_rtdb_get_name(rtdb_name)
      write(*,*) "rtdb_name",rtdb_name
      call util_file_name_resolve(rtdb_name,.false.)
      write(*,*) "resolved rtdb_name",rtdb_name
c
      if(master) then
        inquire(file=rtdb_name,exist=result)
        if(.not.result) then
          call errquit("cannot find rtdb"//
     >                 rtdb_name,
     >                 0,0)
        end if
      end if

      call smd_rtdb_open(rtdb_name,"old")
      call smd_rtdb_get_handle(rtdb)
      write(*,*) "got new rtdb handle",rtdb
      call ga_sync()

      return
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

      subroutine smd_rtdb_get_prefix(prefix)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      character*(*) prefix
c
      integer rtdb
      character*30 pname
      logical result

      pname = "smd_rtdb_get_prefix"
      call smd_rtdb_get_string("smd:prefix",1,prefix,result)
      if(.not.result) call errquit(pname//"no prefix",0,0)


      end 

      subroutine smd_rtdb_get_name1(rtdb_name)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      character*(*) rtdb_name
c
      character*(smd_string_size)  buffer(smd_rtdb_nmax)
      character*(smd_string_size)  tag,sname
      integer rtdb
      character*30 pname
      integer n
      logical result

      pname = "smd_rtdb_get_name1"
      tag = "rtdb"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_data_get_string(sname,"smd:rtdb_name",smd_rtdb_nmax,
     >                          buffer,result)
      if(.not.result) call errquit(pname//"no component",0,0)

      call smd_data_get_int(sname,"smd:rtdb_pointer",1,
     >                         n,result)
      if(.not.result) call errquit(pname//"no component",0,0)

      rtdb_name = buffer(n)

      end 

      subroutine smd_rtdb_push_name(rtdb_name,prefix)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      character*(*) rtdb_name,prefix
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


      call smd_data_get_int(sname,"smd:rtdb_pointer",1,
     >                         n,result)
      if(.not.result) call errquit(pname//"no component",0,0)
c
      n=n+1
      if(n.gt.smd_rtdb_nmax) 
     > call errquit(pname//"out of stack",0,0)
      call smd_data_put_int(sname,"smd:rtdb_pointer",1,n)

      call smd_data_get_string(sname,"smd:rtdb_name",smd_rtdb_nmax,
     >                          buffer,result)
      if(.not.result) call errquit(pname//"no component",0,0)
      buffer(n)=rtdb_name
      call smd_data_put_string(sname,"smd:rtdb_name",smd_rtdb_nmax,
     >                          buffer)

      call smd_data_get_string(sname,"smd:prefix",smd_rtdb_nmax,
     >                          buffer,result)
      if(.not.result) call errquit(pname//"no component",0,0)
      buffer(n)=prefix
      call smd_data_put_string(sname,"smd:prefix",smd_rtdb_nmax,
     >                          buffer)


      end 

      subroutine smd_rtdb_pop_name()
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
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


      call smd_data_get_int(sname,"smd:rtdb_pointer",1,
     >                         n,result)
      if(.not.result) call errquit(pname//"no component",0,0)
c
      n=n-1
      if(n.le.0) 
     > call errquit(pname//"out of stack",0,0)
      call smd_data_put_int(sname,"smd:rtdb_pointer",1,n)

      call smd_data_get_string(sname,"smd:rtdb_name",smd_rtdb_nmax,
     >                          buffer,result)
      if(.not.result) call errquit(pname//"no component",0,0)
      buffer(n+1)=" "
      call smd_data_put_string(sname,"smd:rtdb_name",smd_rtdb_nmax,
     >                          buffer)

      call smd_data_get_string(sname,"smd:prefix",smd_rtdb_nmax,
     >                          buffer,result)
      if(.not.result) call errquit(pname//"no component",0,0)
      buffer(n+1)=" "
      call smd_data_put_string(sname,"smd:prefix",smd_rtdb_nmax,
     >                          buffer)


      end 

      subroutine smd_rtdb_get_prefix1(prefix)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "smd_const_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      character*(*) prefix
c
      character*(smd_string_size)  buffer(smd_rtdb_nmax)
      character*(smd_string_size)  tag,sname
      integer rtdb
      character*30 pname
      integer n
      logical result

      pname = "smd_rtdb_get_prefix1"
      tag = "rtdb"
      call smd_system_get_component(sname,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_data_get_string(sname,"smd:prefix",smd_rtdb_nmax,
     >                          buffer,result)
      if(.not.result) call errquit(pname//"no component",0,0)

      call smd_data_get_int(sname,"smd:rtdb_pointer",1,
     >                         n,result)
      if(.not.result) call errquit(pname//"no component",0,0)

      prefix = buffer(n)

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

