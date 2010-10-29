      program test
C$Id$
      implicit none
c
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer rtdb, ma_handle, ma_index
      integer itest(3), ibuf(3)
      real ftest(4), fbuf(4)
      double precision dtest(5), dbuf(5)
      character *20 cbuf(4), ccbuf(4)
      character *20 name,rtdb_fname
      character *26 date
      logical status
      integer type, nelem, i
c
      data itest/1, 2, 3/
      data ftest/1.0e0, 2.0e0, 3.0e0, 4.0e0/
      data dtest/1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0/
c
      cbuf(1) = 'Have'
      cbuf(2) = 'a'
      cbuf(3) = 'nice'
      cbuf(4) = 'day, Robert!'
c
      call pbeginf()
      if (.not. ma_init(MT_DBL, -1, -1)) stop 99
      call ga_initialize()
      if (.not. rtdb_open('fortran.db', 'empty', rtdb)) stop 1
c
      if (.not. rtdb_put(rtdb, 'Integer array', MT_INT, 3, itest))
     $     stop 2
      if (.not. rtdb_put(rtdb, 'Real array', MT_REAL, 4, ftest))
     $     stop 3
      if (.not. rtdb_put(rtdb, 'Double array', MT_DBL, 5, dtest))
     $     stop 4
      if (.not. rtdb_cput(rtdb, 'Character variable', 4, cbuf))
     $     stop 44
c
      if (.not. rtdb_print(rtdb, .true.)) stop 55
c
      if (.not. rtdb_get(rtdb, 'Integer array', MT_INT, 3, ibuf))
     $     stop 5
      write(6,*) ' Ibuf ', ibuf
      if (.not. rtdb_get(rtdb, 'Real array', MT_REAL, 4, fbuf))
     $     stop 6
      write(6,*) ' Fbuf ', fbuf
      if (.not. rtdb_get(rtdb, 'Double array', MT_DBL, 5, dbuf))
     $     stop 7
      write(6,*) ' Dbuf ', dbuf
      if (.not. rtdb_cget(rtdb, 'Character variable', 4, ccbuf))
     $     stop 8
      write(6,*) ' Cbuf 1 ', ccbuf(1)
      write(6,*) ' Cbuf 2 ', ccbuf(2)
      write(6,*) ' Cbuf 3 ', ccbuf(3)
      write(6,*) ' Cbuf 4 ', ccbuf(4)
c
      status = rtdb_first(rtdb, name)
 10   if (status) then
         if (.not. rtdb_get_info(rtdb, name, type, nelem, date)) stop 9
         write(6,*) name, type, nelem, date
         if (type.eq.mt_int .or. type.eq.mt_dbl) then
            type = -1
            nelem = -1
            if (.not. rtdb_ma_get(rtdb, name, type, nelem, ma_handle))
     $           stop 88
            write(6,*) name, type, nelem, date
            if (.not. ma_get_index(ma_handle, ma_index)) stop 89
            if (type .eq. mt_int) then
               write(6,1) (int_mb(ma_index+i),i=0,nelem-1)
 1             format(1x,10i5)
            else if (type .eq. mt_dbl) then
               write(6,2) (dbl_mb(ma_index+i),i=0,nelem-1)
 2             format(1x,5d16.8)
            else
               stop 91
            endif
            status = ma_free_heap(ma_handle)
         endif
         status = rtdb_next(rtdb, name)
         goto 10
      endif
      if (.not. rtdb_getfname(rtdb, rtdb_fname)) call 
     *     errquit('rtdb_fname call failed',0,0)
          if (.not. rtdb_close(rtdb, 'keep')) stop 10
c
          call testclone(rtdb,rtdb_fname)
c
      call ga_pgroup_set_default(ga_pgroup_get_world())

      call ga_terminate()
      call pend()
c
      end
      subroutine testclone(rtdb,rtdbname)
      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
c
      character*(*) rtdbname           
      character*256 grtdbname
      logical ortdbstate
      double precision dtest(5)
c
      integer ii,ngroups
      integer inode,nnodes,me
      integer maxgroup
      parameter (maxgroup=100)
      integer i,j,proclist(maxgroup),proc_group(0:maxgroup)
      integer mygroup
c
      ngroups=ga_cluster_nnodes()

      if(ga_nodeid().eq.0) write(0,*) '**Testing rtdb_clone'
      ortdbstate=rtdb_parallel(.true.)


      do ii=1,ngroups
         call util_pname0(rtdbname,grtdbname,ii,ngroups)
c     
c     creates rtdb file grtdbname=rtdbname.ii
c     
         if(.not.rtdb_clone(rtdb,grtdbname))
     .        call errquit('oo ',0,0)
      enddo
      me=ga_nodeid()
      inode = ga_cluster_nodeid()
      nnodes = ga_cluster_nprocs(inode)
      if(ga_cluster_nnodes().gt.maxgroup) call 
     g     ga_error('increase maxgroup',-1) 
      do i=0,ga_cluster_nnodes()-1
         do j=0,ga_cluster_nprocs(i)-1
            if(j.gt.maxgroup) call ga_error('increase maxgroup',-1) 
            proclist(j+1)=ga_cluster_procid(i,j)
         enddo
         proc_group(i)=ga_pgroup_create(proclist,ga_cluster_nprocs(i))
      enddo
      call ga_sync()
      do i = 0, ga_cluster_nnodes()-1
         if (i.eq.inode) then
                  mygroup=proc_group(inode)
            write(0,'(I4,A,I4)') me,' joining group',
     ,           proc_group(inode)
            call ga_pgroup_set_default(proc_group(inode))
c     status = ga_create(MT_DBL, n, n, 'a', 1, 1, g_a) 
c     if (.not. status) call ga_error('ga_create failed',-1)
         endif
      enddo
c
c     try  to open rtdbname.igroup
c
      me=ga_nodeid()
      if(me.eq.0) then
         write(0,*) me,' opening rtdb for group ',mygroup
         call util_pname0(rtdbname,grtdbname,mygroup,ngroups)
      endif
      call ga_sync()
      if (.not. rtdb_open(grtdbname, 'old', rtdb)) call errquit
     (     (' rtdb_open old failed ',0,0)
c
c     get the dtest() doubles
c
      if (.not. rtdb_get(rtdb, 'Double array', MT_DBL, 5, dtest))
     $     stop 4
      write(0,123) 'me=',ga_nodeid(),' from group ',mygroup,
     %     '  got dtest(2) =',dtest(2)
c
c     modify dtest(2)
c
      dtest(2)=mygroup*10d0
c
c     rtdb_put for dtest
c
      if (.not. rtdb_put(rtdb, 'Double array', MT_DBL, 5, dtest))
     $   call errquit(' rtdb_put failed ',0,0)
      if (.not. rtdb_get(rtdb, 'Double array', MT_DBL, 5, dtest))
     $     stop 6
      write(0,123) 'me=',ga_nodeid(),' from group ',mygroup,
     ,     '  got dtest(2) =',dtest(2)
      if (.not. rtdb_close(rtdb, 'keep')) call errquit
     (     (' rtdb_close failed ',0,0)
c
c     delete all the rtdb we cloned earlier
c
      if (.not. rtdb_open(grtdbname, 'old', rtdb)) call errquit
     (     (' rtdb_open old failed ',0,0)
         if (.not. rtdb_close(rtdb, 'delete'))call errquit
     (     (' rtdb_close and delete failed ',0,0)
c
 123  format(1X,A,I5,A,I5,A,F10.6)
      return
      end
