C$Id: util_sgroup.F,v 1.3 2004-11-30 01:32:33 edo Exp $
      logical function util_sgstart(rtdb)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
cedo#define GOSMP 1
      integer rtdb
c
      integer nprocs, nsize, midproc, proclist(1000)
      integer i
      character*256 grtdbname
c
c     preliminary rtdb cloning
c
      if (.not. rtdb_getfname(rtdb, rtdb_fname)) call 
     *     errquit('rtdb_fname call failed',0,0)
      if (.not. rtdb_close(rtdb, 'keep')) call errquit(
     *     ' failed to close rtdb for sg ',0,0)
#ifdef GOSMP
      ngroups=ga_cluster_nnodes()
#else
      ngroups=2
#endif
      call util_sgrtdb(rtdb,rtdb_fname,grtdbname)
c
c    Create groups
c   
      write(luout,*) 'Creating groups'
#ifdef GOSMP
      call util_sggosmp()
      util_sgstart=.true.
#else
      nprocs = ga_nnodes()
      midproc = nprocs/2
c
c     hardwired value of 2 subgroups
c
      ngroups=2
      do i = 1, midproc
        proclist(i) = i - 1
      end do
      nsize = midproc
      if (ga_nodeid().eq.0) then
        do i = 1, nsize
          write(luout,*) 'Group 1 ',i,proclist(i)
        end do
      endif
      grouplist(1) = ga_pgroup_create(proclist,nsize)
      do i  = midproc + 1, nprocs
        proclist(i-midproc) = i - 1
      end do
      nsize = nprocs - midproc
      if (ga_nodeid().eq.0) then
        do i = 1, nsize
          write(luout,*) 'Group 2 ',i,proclist(i)
        end do
      endif
      grouplist(2) = ga_pgroup_create(proclist,nsize)
      write(luout,*) 'Groups are created'
      if (ga_nodeid().lt.midproc) then
        call ga_pgroup_set_default(grouplist(1))
c
c     active
c
        util_sgstart=.true.
        mygroup=1
      else
        call ga_pgroup_set_default(grouplist(2))
c
c     inactive
c        
        util_sgstart=.false.
        mygroup=2
      endif
#endif
      sgactive=util_sgstart
      call util_pname0(rtdb_fname,grtdbname,mygroup,ngroups)
      if (.not. rtdb_open(grtdbname, 'old', rtdb)) call errquit
     (     (' rtdb_open old failed ',0,0)
      write(0,*) ga_nodeid(),' of group',mygroup,
     *     'active = ',sgactive, 'rtdbname ',grtdbname
      return
      end
      subroutine util_sgend(rtdb)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
      integer rtdb
c
      if(sgactive) then
         write(0,*) ga_nodeid(),' deleting cloned rtdb '
         if (.not. rtdb_close(rtdb, 'delete'))call errquit
     (        (' sgend: rtdb_close and delete failed ',0,0)
      endif
      write(luout,*) 'end of subgroups call '
      call ga_pgroup_set_default(ga_pgroup_get_world())
c
c     reopen world rtdb
c
      if (.not. rtdb_open(rtdb_fname, 'old', rtdb))
     $     call errquit('sgend: rtdb_open old failed ', 0, 0)
      
      return
      end
      subroutine util_sgrtdb(rtdb,rtdbname,grtdbname)
      implicit none
      integer rtdb
      character*(*) rtdbname,grtdbname
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
      integer ii
      logical ortdbstate,ortdbstate2
c
      write(luout,*) 'util_sgrtdb start'
      do ii=1,ngroups
         call util_pname0(rtdbname,grtdbname,ii,ngroups)
c     
c     creates rtdb file grtdbname=rtdbname.ii
c     
         ortdbstate=rtdb_parallel(.true.)
         write(0,*) ' rtdb clone',rtdb,grtdbname
         if(.not.rtdb_clone(rtdb,grtdbname))
     .        call errquit('oo ',0,0)
         ortdbstate2=rtdb_parallel(ortdbstate)
      enddo
      call ga_sync()
      return
      end
      subroutine util_sggosmp()
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "util_sgroup.fh"
      integer me,inode,nnodes
      integer i,j,proclist(maxgroups)
c
      ngroups=ga_cluster_nnodes()

      me=ga_nodeid()
      inode = ga_cluster_nodeid()
      nnodes = ga_cluster_nprocs(inode)
      if(ga_cluster_nnodes().gt.maxgroups) call 
     g     errquit('sggosmp: increase maxgroups',0,0) 
      do i=0,ga_cluster_nnodes()-1
         do j=0,ga_cluster_nprocs(i)-1
            proclist(j+1)=ga_cluster_procid(i,j)
         enddo
         grouplist(i+1)=ga_pgroup_create(proclist,ga_cluster_nprocs(i))
      enddo
      call ga_sync()
      do i = 0, ga_cluster_nnodes()-1
         if (i.eq.inode) then
                  mygroup=grouplist(inode+1)
            write(0,'(I4,A,I4)') me,' joining group',
     ,          grouplist(inode+1)
            call ga_pgroup_set_default(mygroup)
         endif
      enddo
      return
      end
