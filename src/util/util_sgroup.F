C$Id: util_sgroup.F,v 1.7 2004-12-01 18:27:38 edo Exp $
      logical function util_sgstart(rtdb)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
#include "mafdecls.fh"
      integer rtdb
c
      integer nprocs, nsize, midproc, proclist(1000)
      integer i,errout
      character*256 grtdbname
      logical parallel, gosmp, lstatus
c
      if (.not.rtdb_get(rtdb, 'rtdbio_parallel', mt_log, 1,
     &   parallel))parallel=.false.
      if (.not.rtdb_get(rtdb, 'gosmp', mt_log, 1,
     &   gosmp))gosmp=.false.
c
c     preliminary rtdb cloning
c
      if (.not. rtdb_getfname(rtdb, rtdb_fname)) call 
     *     errquit('rtdb_fname call failed',0,0)
      if (.not. rtdb_close(rtdb, 'keep')) call errquit(
     *     ' failed to close rtdb for sg ',0,0)
c
c     predefine ngroups for rtdb cloning
c
      if(gosmp) then
!         ngroups=ga_cluster_nnodes()
      else
         ngroups=2
      endif
c
c    Create groups
c   
      write(luout,*) 'Creating groups'
      if(gosmp) then
         call util_sggosmp()
         call util_sgrtdb(rtdb,rtdb_fname,grtdbname,parallel)
         util_sgstart=.true.
      else
         call util_sgrtdb(rtdb,rtdb_fname,grtdbname,.true.)
         call util_sggoidle(util_sgstart)
      endif
      sgactive=util_sgstart
      if(parallel) then
      call util_pname0(rtdb_fname,grtdbname,mygroup,ngroups)
      endif
      if (.not. rtdb_open(grtdbname, 'old', rtdb)) call errquit
     (     (' rtdb_open old failed ',0,0)
      write(0,*) ga_nodeid(),' of group',mygroup,
     *     'active = ',sgactive, 'rtdbname ',grtdbname
c
c     nasty hack stdoutput redirect to fort.90+group for 
c     groupid gt 0
c
      luout_0=luout
      if(mygroup.ne.1.and.gosmp) then
         luout=90+mygroup
         write(0,*) mygroup,'new value for luout',luout
         if(ga_nodeid().eq.0) then
            inquire(luout, exist=lstatus,iostat=errout)
            if(lstatus) then
               open(luout,  status='old',iostat=errout)
               if(errout.eq.0)close(luout, status='delete')
            else
               write(0,*) ' fort.9x not there'
            endif
            open(luout, status='new')
         endif
         call ga_sync()
      endif
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
      luout=luout_0
      write(luout,*) 'end of subgroups call '
      call ga_pgroup_set_default(ga_pgroup_get_world())
c
c     reopen world rtdb
c
      if (.not. rtdb_open(rtdb_fname, 'old', rtdb))
     $     call errquit('sgend: rtdb_open old failed ', 0, 0)
      
      return
      end
      subroutine util_sgrtdb(rtdb,rtdbname,grtdbname,
     ,     parallel)
      implicit none
      integer rtdb
      character*(*) rtdbname,grtdbname
      logical parallel
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util_sgroup.fh"
      integer ii
      logical ortdbstate,ortdbstate2,status
      character*256 scratchrtdb,cmd,cmd1
      character*10 pname
      integer util_system
      external util_system
c
      write(luout,*) 'util_sgrtdb start'
      if(parallel) then
         do ii=1,ngroups
            call util_pname0(rtdbname,grtdbname,ii,ngroups)
c     
c     creates rtdb file grtdbname=rtdbname.ii
c     
            ortdbstate=rtdb_parallel(.true.)
            if(.not.rtdb_clone(rtdb,grtdbname))
     .           call errquit('oo ',0,0)
            ortdbstate2=rtdb_parallel(ortdbstate)
         enddo
         call ga_sync()
      else
c
c     do a copy of rtdb to scratch_dir for each SMP node
c
         write(pname,333) mygroup
 333       format('sg_db.',i4.4)

            call util_file_name(pname, .true., .false., 
     $           scratchrtdb)
c
c     assuming we are already in subgroups ...
c     am I proc 0 of a subgroup ?
c
         if(ga_nodeid().eq.0) then
c
c     test if world rtdb is around
c
            inquire(file=rtdbname(1:inp_strlen(rtdbname)),
     *       exist=status)
            if(.not.status) then
c
c     most likely cause of error: permdir not available on 
c     this node
c
               write(0,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               write(0,*) ' ! please change permanent_dir to a   !'
               write(0,*) ' ! filesystem available to all the    !'
               write(0,*) ' ! processors (e.g. NFS, GPFS, Lustre)!'
               write(0,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
               call errquit(
     *              ' util_sgrtdb: failed to open rtdb',0,0)
            endif
            cmd = 'cp -f '
            write(cmd1,'(a,'' '',a,'' '',a)')
     *           cmd(1:inp_strlen(cmd)),
     *           rtdbname(1:inp_strlen(rtdbname)),
     *           scratchrtdb(1:inp_strlen(scratchrtdb))
            cmd = cmd1
            if (util_system (cmd) .eq. 0) then
               Write (0,*) cmd, ' OK'
            else
               Write (0,*) cmd, ' failed '
!               call errquit('util_sgrtdb: problem with system call',
!     *              0, 0)
            endif
         endif
            grtdbname=scratchrtdb
      endif
      return
      end
      subroutine util_sggosmp()
      implicit none
c
c     creates one subgroup for each smp
c
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
      subroutine util_sggoidle(active)
      implicit none
c
c     creates two subgroups
c     one active
c     second idle
c
#include "stdio.fh"
#include "global.fh"
#include "util_sgroup.fh"
      logical active
c
      integer nproc,midproc,nprocs
      integer i,proclist(maxgroups)
      integer nsize
c      
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
        active=.true.
        mygroup=1
      else
        call ga_pgroup_set_default(grouplist(2))
c
c     inactive
c        
        active=.false.
        mygroup=2
      endif
      end
