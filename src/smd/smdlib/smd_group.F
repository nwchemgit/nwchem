      block data initial_smd_group_data
      implicit none
#include "smd_group_data.fh"
c     
      data  smd_group_luout /6/
      data  smd_ngroup /-1/
      data  smd_group_nproc /-1/
      data  smd_group_id /-1/
      data  smd_group_handle /-1/
      data smd_group_rtdb_reuse /.false./

      end
      subroutine smd_group_init(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
c     
      integer rtdb
c
      integer ngroups
      integer i1

      character*255 buffer, buffer1
      if (.not.rtdb_get(rtdb, 'subgroups_number', mt_int, 1,
     &   ngroups)) 
     $  call errquit('define number of groups', 0, RTDB_ERR)
      call ga_sync()

      if(.not.rtdb_getfname(rtdb,smd_group_rtdb_world))
     $  call errquit('rtdb_getfname', 0, RTDB_ERR)

      call ga_sync()
      if(.not.rtdb_close(rtdb, 'keep'))
     + call errquit('Failed to close group rtdb',0, GEOM_ERR)

      call smd_group_create_simple(ngroups)

c      call util_file_prefix_get(buffer)
c      i1 = inp_strlen(buffer)
c      write(buffer1,'(A,I4.4)')
c     >    buffer(1:i1),smd_group_id
c      call util_file_prefix_set(buffer1)
c
c      call util_file_prefix_get(smd_group_prefix_world)
c      i1 = inp_strlen(smd_group_prefix_world)
c      write(smd_group_prefix,'(A,I4.4)')
c     >    smd_group_prefix_world(1:i1),smd_group_id
c      call util_file_prefix_set(smd_group_prefix)


      return
      end
      subroutine smd_group_set_io()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
#include "stdio.fh"
c     
      integer fn,i1
      logical result,ogroupio
c
      character*255 buffer,buffer1

#ifdef USE_SUBGROUPS      
      call smd_rtdb_get_log("smd:groupio",1,ogroupio,result)
      if(.not.result) return
      if(.not.ogroupio) return
      call util_file_prefix_get(buffer)
      i1 = inp_strlen(buffer)
      write(buffer1,'(A,I4.4,I4.4,".out")')
     >    buffer(1:i1),smd_group_id,ga_nodeid()
      
      if(.not.util_get_io_unit(fn))
     + call errquit('no units available',0, 0)

      close(luout)
      luout = fn
      open(unit=luout,status="unknown",form="formatted",file=buffer1)

#endif
      return
      end
      subroutine smd_group_set_io_custom(im,stub)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
#include "stdio.fh"
      integer im
      character*(*) stub
c     
      integer fn,i1,i2
      integer iw
c
      character*255 buffer,buffer1
      character*10 fs
      logical ogroupio,result

#ifdef USE_SUBGROUPS      
      call smd_rtdb_get_log("smd:groupio",1,ogroupio,result)
      if(.not.result) return
      if(.not.ogroupio) return
 
      write(fs,'(I10)') im
      call inp_adjustl(fs)
      i1 = inp_strlen(stub)
      i2 = inp_strlen(fs)
      write(buffer1,'(A,A,".out")')
     >    stub(1:i1),fs(1:i2)
      
      if(.not.util_get_io_unit(fn))
     + call errquit('no units available',0, 0)

c
c     leave unit 6 alone
c     ------------------
c      if(luout.ne.6) then
c        close(luout)
c      end if
      luout = fn
      open(unit=luout,status="unknown",position="append",
     >       form="formatted",file=buffer1)
      write(luout,*) "Calculation from group",smd_group_id


#endif
      return
      end

      subroutine smd_group_unset_io()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
#include "stdio.fh"
c     
      integer fn,i1
c
      character*255 buffer,buffer1
      logical result,ogroupio

      call smd_rtdb_get_log("smd:groupio",1,ogroupio,result)
      if(.not.result) return
      if(.not.ogroupio) return

      close(luout)

      return
      end
      subroutine smd_group_prefix_set(im)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
c     
      integer im
      integer rtdb
c
      character*255 buffer,buffer1
      integer i1

      call util_file_prefix_get(buffer)
      i1 = inp_strlen(buffer)
      write(buffer1,'(A,I4.4)')
     >    buffer(1:i1),im
      call util_file_prefix_set(buffer1)

      if(ga_nodeid().eq.0)
     > call util_flush(6)
      call ga_sync(ga_pgroup_get_world())
      if(smd_group_id.eq.1.and.ga_nodeid().eq.0) then
        write(*,*) "GROUP INFORMATION"
        call util_flush(6)
      end if
      if(ga_nodeid().eq.0) then
       write(*,*) "Group id, nprocs",smd_group_id,smd_group_nproc
       call util_flush(6)
      end if
      call ga_sync(ga_pgroup_get_world())

      return
      end
      subroutine smd_group_end()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
c     
      integer rtdb
c
      integer ngroups

      character*36 buffer

      call ga_pgroup_set_default(ga_pgroup_get_world())
      if(.not. ga_pgroup_destroy(smd_group_handle))
     + call errquit('Failed to destroy group',0, 0)
      call ga_sync()

      smd_ngroup =-1
      smd_group_nproc =-1
      smd_group_id =-1
      smd_group_handle =-1

      return
      end

      subroutine smd_group_rtdb_clone1()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
c     
      character*72 buffer
      integer i1,i2
      logical result
      buffer = " "
c      i1 = 0
c10    if (inp_strtok(smd_group_rtdb_world, '/',i1,i2)) then
c      buffer = smd_group_rtdb_world(i1:i2)
c      goto 10
c      endif
c
c      write(smd_group_rtdb,'(a,I4.4)')
c     >          buffer(1:inp_strlen(buffer)),
c     >          smd_group_id

      write(smd_group_rtdb,'(I4.4,A3)')
     >          smd_group_id,".db"

      buffer = smd_group_rtdb
      call util_file_name(buffer,.true.,.false.,smd_group_rtdb)
      if(ga_nodeid().eq.0) then
        inquire(file=smd_group_rtdb_world,exist=result)
        if(.not.result) then
          call errquit("cannot find world rtdb"//
     >                 smd_group_rtdb_world,
     >                 0,0)
        end if
        call util_file_copy(smd_group_rtdb_world,
     >                      smd_group_rtdb)
        inquire(file=smd_group_rtdb,exist=result)
        if(.not.result) then
          call errquit("cannot find group rtdb"//
     >                 smd_group_rtdb,
     >                 0,0)
        end if
 
      end if
      return
      end

      subroutine smd_group_rtdb_clone2(rtdb,im)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
#include "rtdb.fh"
      integer rtdb
      integer im
c     
      character*255 buffer,buffer1
      integer i1,i2
      logical result
c
      call ga_sync()
c
      buffer = ""
      write(buffer,'(I4.4,A3)')
     >          im,".db"

      call util_file_prefix(buffer,smd_group_rtdb)
      if(ga_nodeid().eq.0) then
        if(smd_group_rtdb_reuse) then
        buffer = smd_group_rtdb
        call util_file_name_resolve(buffer,.false.)
        inquire(file=buffer,exist=result)
        else
        result = .false.
        end if
        if(.not.result) then
          buffer = smd_group_rtdb_world
        end if
        buffer1 = smd_group_rtdb
        call util_file_name_resolve(buffer1,.true.)
        write(*,*) "copying rtdb from",buffer(1:inp_strlen(buffer))
        call util_file_copy(buffer,buffer1)
      end if

      if(.not.rtdb_open(buffer1,"old",rtdb))
     + call errquit('Failed to open group rtdb',0, GEOM_ERR)
      
      smd_group_rtdb_reuse = .true.

      return
      end

      subroutine smd_group_rtdb_clone(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
#include "rtdb.fh"
      integer rtdb
c     
      character*72 buffer
      integer i1,i2
      logical result
c
      call ga_sync()
c
      write(smd_group_rtdb,'(I4.4,A3)')
     >          smd_group_id,".db"

      buffer = smd_group_rtdb
      call util_file_name(buffer,.true.,.false.,smd_group_rtdb)
      if(ga_nodeid().eq.0) then
        inquire(file=smd_group_rtdb_world,exist=result)
        if(.not.result) then
          call errquit("cannot find world rtdb"//
     >                 smd_group_rtdb_world,
     >                 0,0)
        end if
        call util_file_copy(smd_group_rtdb_world,
     >                      smd_group_rtdb)
        inquire(file=smd_group_rtdb,exist=result)
        if(.not.result) then
          call errquit("cannot find group rtdb"//
     >                 smd_group_rtdb,
     >                 0,0)
        end if
 
      end if

      if(.not.rtdb_open(smd_group_rtdb,"old",rtdb))
     + call errquit('Failed to open group rtdb',0, GEOM_ERR)

      return
      end

      subroutine smd_group_rtdb_release(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
#include "rtdb.fh"
      integer rtdb
c     
      character*72 buffer
      integer i1,i2
      logical result
c
      if(.not.rtdb_close(rtdb, 'delete'))
     + call errquit('Failed to close group rtdb',0, GEOM_ERR)
c
      smd_group_rtdb = " "


      return
      end

      subroutine smd_group_rtdb_release1(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "smd_group_data.fh"
#include "rtdb.fh"
      integer rtdb
c     
      character*255 buffer,buffer1
      integer i1,i2
      logical result

      if(.not.rtdb_getfname(rtdb,buffer))
     $  call errquit('rtdb_getfname', 0, RTDB_ERR)
c
c
      if(ga_nodeid().eq.0) then
        buffer1 = smd_group_rtdb
        call util_file_name_resolve(buffer1,.false.)
        call util_file_copy(buffer,buffer1)
      end if

c
      if(.not.rtdb_close(rtdb, 'delete'))
     + call errquit('Failed to close group rtdb',0, GEOM_ERR)
c
      smd_group_rtdb = " "


      return
      end

      subroutine smd_group_inquire(rtdb,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      logical util_sgstart
      external util_sgstart
c
      integer rtdb
      logical result
c
      integer ngroup

      if (.not.rtdb_get(rtdb, 'subgroup_number', mt_int, 1,
     &   ngroup)) ngroup = 1

      result = ngroup.gt.1

      end

      subroutine smd_group_create_simple0(group_want)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_group_data.fh"
c
c     Create subgroup of a constant size
c
      integer group_want
      integer nproc,myproc,nchunkq,nremainq,i,j,n
      integer num_proc
      integer i_l,nl,h_l, plo, phi

      character*30 pname
c
      pname = "smd_group_create_simple"
c 
      smd_ngroup = group_want
c
      nproc = ga_nnodes()
      myproc = ga_nodeid()

      nchunkq=nproc/smd_ngroup
      nremainq=mod(nproc,smd_ngroup)

      nl = nchunkq+nremainq
      if(.not.ma_push_get(mt_int,nl,'proclist',h_l,i_l))
     + call errquit(pname//'Failed to allocate memory for i_l',
     + nl, MA_ERR)

      call ifill(nl,-1,int_mb(i_l),1)


      n = 0
      plo = 0
      do i=1,smd_ngroup
         num_proc = nchunkq
         if(i .le. nremainq) num_proc = num_proc + 1
         phi = plo + num_proc
         do j=1,num_proc
            int_mb(i_l+j-1) = plo+j
         enddo
         if ( myproc.ge.plo .and. myproc.lt.num_proc) then
            smd_group_handle = ga_pgroup_create(int_mb(i_l), num_proc)
            smd_group_nproc = num_proc
         endif
      enddo

c set the group

      if(.not.ma_pop_stack(h_l))
     & call errquit(pname//'Failed to deallocate stack h_itmp',0,
     &       MA_ERR)

      call ga_pgroup_set_default(smd_group_handle)

      return
      end

      subroutine smd_group_create_simple(group_want)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_group_data.fh"
c
c     Create subgroup of a constant size
c
      integer group_want
      integer nproc,myproc,nchunkq,nremainq,i,j,n
      integer num_proc
      integer i_l,nl,h_l

      character*30 pname
c
      pname = "smd_group_create_simple"
c 
      smd_ngroup = group_want
c
      nproc = ga_nnodes()
      myproc = ga_nodeid()

      nchunkq=nproc/smd_ngroup
      nremainq=mod(nproc,smd_ngroup)

      nl = nchunkq+nremainq
      if(.not.ma_push_get(mt_int,nl,'proclist',h_l,i_l))
     + call errquit(pname//'Failed to allocate memory for i_l',
     + nl, MA_ERR)

      call ifill(nl,-1,int_mb(i_l),1)

      n = 0
      do i=1,smd_ngroup
         num_proc = nchunkq
         if(i .le. nremainq) num_proc = num_proc + 1
         do j=1,num_proc
            int_mb(i_l+j-1)=n
            if (n .eq. myproc) smd_group_id = i
            n=n+1
         enddo
         if (smd_group_id .eq. i) then
            smd_group_handle = ga_pgroup_create(int_mb(i_l), num_proc)
            smd_group_nproc = num_proc
         endif
      enddo
 
c      write(*,*) "processor ",myproc, "blongs to ",smd_group_id 
      call util_flush(6)
      call ga_sync()

      if(.not.ma_pop_stack(h_l))
     & call errquit(pname//'Failed to deallocate stack h_itmp',0,
     &       MA_ERR)

c set the group
      smd_group_luout = 170+smd_group_id
      call ga_pgroup_set_default(smd_group_handle)
c
c     this call has been now moved to task initialization part
c     --------------------------------------------------------
c      call smd_group_set_io()
      return
      end

      subroutine smd_group_ngroup_get(ngroup)
      implicit none
#include "smd_group_data.fh"
      integer  ngroup
      ngroup = smd_ngroup
      return
      end

      subroutine smd_group_id_get(group_id)
      implicit none
#include "smd_group_data.fh"
      integer  group_id
      group_id = smd_group_id
      return
      end

      subroutine smd_group_luout_get(group_luout)
      implicit none
#include "smd_group_data.fh"
      integer  group_luout
      group_luout = smd_group_luout
      return
      end

      subroutine smd_group_nproc_get(group_nproc)
      implicit none
#include "smd_group_data.fh"
      integer  group_nproc
      group_nproc = smd_group_nproc
      return
      end

      subroutine smd_group_message(buffer)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "smd_group_data.fh"
      character*(*) buffer
c
      integer un

      if(ga_nodeid().eq.0) then
      un = 160+smd_group_id
      write(un,*) buffer,ga_nodeid()
      call util_flush(un) 
      end if
      call ga_sync()
      return
      end

c $Id$
