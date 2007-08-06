C$Id: util_sgroup.F,v 1.34 2007-08-06 21:00:22 d3p852 Exp $

C
      subroutine util_sgstart2(groups_want,rtdb)
C     This is a python interface routine
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb, groups_want
      logical util_sgstart
      external util_sgstart
      logical ignore

      ignore=rtdb_put(rtdb, 'subgroups_number', mt_int, 1, groups_want)
      ignore=util_sgstart(rtdb)

      end 


      logical function util_sgstart(rtdb)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
      integer rtdb
c
      integer nprocs, nsize, midproc
      integer i,errout
      character*256 grtdbname
      logical lstatus

C strings for fname definition
      character*256 permdir
      character*256 fprefix
      integer myproc,nproc
c
c Read number of subgroups from input
c
      if (ga_pgroup_get_default() .eq. ga_pgroup_get_world()) then
        depth = 0;
        mygroup(depth) = 1
        ngroups(depth) = 1
        my_ga_grp(depth) = ga_pgroup_get_world()
      end if
      depth = depth + 1
      if (depth .gt. maxdeep) 
     *     call errquit('depth is greater than maxdeep',0,0)

      myproc=ga_nodeid()
      nproc=ga_nnodes()
c this is added here to print out sub groups info. !DW
      sgprint = util_print('sgroup info', print_debug)

      if (.not.rtdb_get(rtdb, 'subgroups_number', mt_int, 1,
     &   ngroups(depth)))ngroups(depth)=1

c
c     preliminary rtdb cloning
c
C  Trying to fix the 36 char cutoff for fname
      if (.not. rtdb_getfname(rtdb, rtdb_fname(depth))) call 
     *     errquit('rtdb_fname call failed',0,0)
C  Added to overcome 36 character limit
      write(permdir,'(256(a))') (' ', i=1,256)
      write(fprefix,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'permanent_dir',1,permdir))
     &     call errquit('rtdb get permanent_dir failed',0,0)
      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix))
     &     call errquit('rtdb get file_prefix failed',0,0)
      write(rtdb_fname(depth),'(a,a,a,a)')
     &  permdir(15:inp_strlen(permdir)),
     &  '/',
     &  fprefix(1:inp_strlen(fprefix)),
     &  '.db'
C  End of Added fix
      if (.not. rtdb_close(rtdb, 'keep')) call errquit(
     *     ' failed to close rtdb for sg ',0,0)
c
c    Create groups
c   
      if(myproc.eq.0) write(luout,*) 'Creating groups'

c initialize sub groups using one general method

      call util_sginit()
C     call util_sginit_smp()
C     call util_sginit_irreg(cpusperg_array)

      call util_sgrtdb(rtdb,rtdb_fname(depth),grtdbname)
      util_sgstart=.true.
      if (.not. rtdb_open(grtdbname, 'old', rtdb)) call errquit
     *     (' rtdb_open old failed ',0,0)

      if(myproc.eq.0.and.sgprint) then
        write(6,*) ga_nodeid(),' of group',mygroup(depth),'rtdbname '
     *           ,grtdbname, ' rtdb=',rtdb
      endif
c
c     Create movecs for each subgroup
c
      call util_sgmovecs(rtdb)
c
c
      luout_0(depth)=luout
      if(myproc.eq.0.and.sgprint) then
        write(6,*) 'everyone should have an open file at this point'
      endif
      return
      end

      subroutine util_sgend(rtdb)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
      integer rtdb,myproc
      logical ignore
      integer mypgroup_old
c
      myproc=ga_nodeid()
      if(myproc.eq.0.and.sgprint) then
        write(0,*) ' deleting cloned rtdb '
      endif
      if (.not. rtdb_close(rtdb, 'delete'))call errquit
     $        (' sgend: rtdb_close and delete failed ',0,0)
      if(mygroup(depth).ne.1) then
        close(luout, status='keep')
      endif
      luout=luout_0(depth)
      if(myproc.eq.0) then
        write(luout,*) 'end of subgroups call '
      endif
      mypgroup_old = ga_pgroup_get_default()
      call ga_pgroup_set_default(my_ga_grp(depth-1))
      ignore = ga_pgroup_destroy(mypgroup_old)
c
c     reopen world rtdb
c
      if (.not. rtdb_open(rtdb_fname(depth), 'old', rtdb))
     $     call errquit('sgend: rtdb_open old failed ', 0, 0)
      
      depth = depth - 1

      return
      end


      subroutine util_sgrtdb(rtdb,rtdbname,grtdbname)
      implicit none
      integer rtdb
      character*(*) rtdbname,grtdbname
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util_sgroup.fh"
      integer ii
      logical ortdbstate,ortdbstate2,status
      character*256 scratchrtdb
      character*10 pname
      integer util_system
      external util_system
      external util_file_copy
      integer myproc
c
      myproc=ga_nodeid()

      if(myproc.eq.0.and.sgprint) then
        write(luout,*) 'util_sgrtdb start'
      endif
c
c     do a copy of rtdb to scratch_dir for each SMP node
c
      write(pname,333) mygroup(depth)
 333     format('sg_db.',i4.4)

      call util_file_name(pname, .true., .false., scratchrtdb)
c
c     assuming we are already in subgroups ...
c     am I proc 0 of a subgroup ?
c
      if(myproc.eq.0.and.sgprint) then
         write(0,*)'Entering util_sgrtdb if(node 0)',ga_nodeid(),
     &   ga_pgroup_get_default()
      endif
      call flush(0)
      if(ga_nodeid().eq.0) then
c
c     test if world rtdb is around
c
        if(myproc.eq.0.and.sgprint) then
           write(6,*)'world rtdbname is ',
     &       rtdbname(1:inp_strlen(rtdbname))
           write(6,*)'group rtdbname is ',
     &       scratchrtdb(1:inp_strlen(scratchrtdb))
        endif
        inquire(file=rtdbname(1:inp_strlen(rtdbname)),exist=status)
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
           call errquit(' util_sgrtdb: failed to open rtdb',0,0)
        endif
        call util_file_copy(rtdbname(1:inp_strlen(rtdbname)),
     &        scratchrtdb(1:inp_strlen(scratchrtdb)))
C  Trying to remove system call (Crashes MVAPICH)
        if(myproc.eq.0.and.sgprint) then
            write (0,*)'rtdb copied',ga_nodeid(),ga_pgroup_get_default()
        endif
        inquire(file=scratchrtdb(1:inp_strlen(scratchrtdb)),
     &          exist=status)
        if (status) then
          if(myproc.eq.0.and.sgprint) then
            Write (0,*)'copy OK'
          endif
        else
          Write (0,*) 'copy failed'
          call errquit('util_sgrtdb: problem with system call', 0, 0)
        endif
      endif
      grtdbname=scratchrtdb
      if(myproc.eq.0) print*,'rtdb copied'
      return
      end

      subroutine util_sginit_smp()
      implicit none
c
c     creates one subgroup for each smp
c
#include "global.fh"
#include "util_sgroup.fh"
      integer inode
      integer i,j,proclist(maxcpus)
      integer grouplist(maxgroups)

c setup

      inode = ga_cluster_nodeid()
      ngroups(depth) = ga_cluster_nnodes()

      if(ngroups(depth) .gt.maxgroups)
     *      call errquit('increase maxgroups',0,0) 

c make each smp node into a group

      do i=0,ngroups(depth)-1
         do j=0,ga_cluster_nprocs(i)-1
           if(j+1 .gt. maxcpus)
     *            call errquit('increase maxcpus',0,0)
            proclist(j+1)=ga_cluster_procid(i,j)
         enddo
         grouplist(i+1)=ga_pgroup_create(proclist,ga_cluster_nprocs(i))
      enddo

C Set the group

      mygroup(depth)=inode+1
      my_ga_grp(depth) = grouplist(mygroup(depth))
      write(0,'(I4,A,I4)')ga_nodeid(),' joining group', my_ga_grp(depth)
      call ga_pgroup_set_default(my_ga_grp(depth))

      return
      end

      subroutine util_sginit_irreg(cpusperg_array)
      implicit none
c
c     creates custom subgroups using cpusperg_array()
c
#include "global.fh"
#include "util_sgroup.fh"
      integer inode
      integer i,j,proclist(maxcpus), myprocid,proc_counter,group_counter
      integer grouplist(maxgroups)
      integer cpusperg_array(maxgroups)
c
      inode = ga_cluster_nodeid()

C This code loops over processes and divides them into groups based 
C upon cpusperg_array() and stores this information in proclist().
C At the same time, it figures out which group the local process is in
C and stores this information in mygroup(depth).
      proc_counter=1
      group_counter=1
      do i=0,ga_cluster_nnodes()-1
         do j=0,ga_cluster_nprocs(i)-1
           if(proc_counter .gt. maxcpus)
     *            call errquit('increase maxcpus',0,0)
           proclist(proc_counter)=ga_cluster_procid(i,j)
           myprocid=ga_cluster_procid(inode,j)
           if (proclist(proc_counter).eq.myprocid) 
     $          mygroup(depth)=group_counter
           if (cpusperg_array(group_counter).eq.proc_counter) then
             if(group_counter.gt.maxgroups) 
     *                 call errquit('increase maxgroups',0,0) 
             grouplist(i)=ga_pgroup_create(proclist,proc_counter)
             group_counter=group_counter+1
             proc_counter=0
           end if
           proc_counter=proc_counter+1
         enddo
      enddo

      ngroups(depth) = group_counter - 1
      my_ga_grp(depth) = grouplist(mygroup(depth))
      call ga_pgroup_set_default(my_ga_grp(depth))

      write (6,*) 'me=',ga_nodeid(),' is the part of subgroup=',
     *            mygroup(depth)

      return
      end

      subroutine util_sginit()
      implicit none
#include "global.fh"
#include "util_sgroup.fh"
c
c     Create subgroups of a constant size
c
      integer nproc,myproc,nchunkq,nremainq,i,j,n,meafter
      integer num_proc,proclist(maxcpus), grouplist(maxgroups)

c  on either smp or cluster of smp, ga_nnodes is consistent value, 
c  so ga_nnodes as starting base
      nproc = ga_nnodes()
      myproc = ga_nodeid()

c split number of procs into groups, starting group index:0
      nchunkq=nproc/ngroups(depth)
      nremainq=mod(nproc,ngroups(depth))
      if (nremainq .eq. 0) then
        if(nchunkq .gt. maxcpus)
     *      call errquit('increase maxcpus',0,0)
      else
        if(nchunkq+1 .gt. maxcpus)
     *      call errquit('increase maxcpus',0,0)
      end if
      if(ngroups(depth).gt.maxgroups)
     *      call errquit('increase maxgroups',0,0)

c for irregular distribution, make sure the lower groups have more procs
c assign each proc to individual group, then create individual group
      n = 0
      do i=1,ngroups(depth)
         num_proc = nchunkq
         if(i .le. nremainq) num_proc = num_proc + 1
         do j=1,num_proc
            proclist(j)=n
            if (n .eq. myproc) mygroup(depth) = i
            n=n+1
         enddo
         grouplist(i)=ga_pgroup_create(proclist, num_proc)
      enddo

c set default individual group
      my_ga_grp(depth) = grouplist(mygroup(depth))
      call ga_pgroup_set_default(my_ga_grp(depth))

      meafter=ga_nodeid()
      
      write(*,"(2x,'proc.',I4,' now is proc.',I4,' in group',I4)")
     $myproc,meafter,mygroup(depth)

      return
      end
c
      subroutine util_sgmovecs(rtdb)
      implicit none
      integer rtdb
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util_sgroup.fh"
      character*255 movecs_in,sg_movecs_out,tag,theory
      logical movecs_out_l
      integer myproc
c
c Tested only for the case when movecs are created
c

c
c Read movecs filename with path
c
      myproc=ga_nodeid()

       if (.not.(rtdb_cget(rtdb,'scf:input vectors',1,movecs_in).or.
     $           rtdb_cget(rtdb,'dft:input vectors',1, movecs_in)))
     $     movecs_in = 'atomic'
       if (movecs_in.eq.'atomic')
     $     call util_file_name('movecs',.true.,.false.,movecs_in)
c
c Create movecs.mygroup(depth)
c
      call util_pname0(movecs_in,sg_movecs_out,mygroup(depth),
     *     ngroups(depth))
  
      if(myproc.eq.0.and.sgprint) then
        write(6,*) 'mygroup=',mygroup(depth),
     $             ' movecs_in',movecs_in,
     $             ' sg_movecs_out=',sg_movecs_out
      endif
c
c Set output movecs
c
      if (.not.(rtdb_cget(rtdb,'task:theory',1,theory)))
     *     call errquit('rtdb_get task:theory failed',0,0)
      if (theory.eq.'dft') then
       tag='dft:output vectors'
      else
       tag='scf:output vectors'
      endif
c
      if (.not. rtdb_cput(rtdb,tag,1,sg_movecs_out)) then
       write (6,*) 'error in rtdb_cput'
      endif
c
      return
      end

c     Hokey functions to keep others away from my common blocks

      integer function util_sgroup_mygroup()
#include "util_sgroup.fh"
      util_sgroup_mygroup = mygroup(depth)
      end
      integer function util_sgroup_numgroups()
#include "util_sgroup.fh"
      util_sgroup_numgroups = ngroups(depth)
      end

