C$Id: util_sgroup.F,v 1.20 2007-08-03 22:47:11 d3p852 Exp $

C
C TODO:  Need to be able to subdivide existing groups
C TODO:  Need to reset the number of groups in sgend
C TODO:  Need to clean up  references to world
C        Right now, if you have python generate groups and then have
C          a method that uses groups, you are hosed, and very strange
C          things will happen (hopefully you with crash).
C WARNING:  mygroup is NOT the GA group, but rather a 1 to N
C           numbering scheme.
C

      subroutine util_sgstart2(groups_want,rtdb)
C     This is a python interface routine
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
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
      integer nprocs, nsize, midproc, proclist(1000)
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
c util_sginit is rewritten to allow creating non uniform size of procs 
c in subgroups on both SMP and cluster of SMP machines !DW

C TODO: At some point this would be nice to be able to do
      if (ga_pgroup_get_default() .ne. ga_pgroup_get_world()) then
         write(0,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write(0,*) ' ! Cannot create groups when already  !'
         write(0,*) ' ! broken into groups                 !'
         write(0,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         call errquit(' util_sgstart: failed to make groups',0,0)
      end if

      myproc=ga_nodeid()
      nproc=ga_nnodes()
c this is added here to print out sub groups info. !DW
      sgprint = util_print('sgroup info', print_debug)

      if (.not.rtdb_get(rtdb, 'subgroups_number', mt_int, 1,
     &   ngroups))ngroups=1

c
c     preliminary rtdb cloning
c
C  Trying to fix the 36 char cutoff for fname
      if (.not. rtdb_getfname(rtdb, rtdb_fname)) call 
     *     errquit('rtdb_fname call failed',0,0)
C  Added to overcome 36 character limit
      write(permdir,'(256(a))') (' ', i=1,256)
      write(fprefix,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'permanent_dir',1,permdir))
     &     call errquit('rtdb get permanent_dir failed',0,0)
      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix))
     &     call errquit('rtdb get file_prefix failed',0,0)
      write(rtdb_fname,'(a,a,a,a)')
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
      call util_sgrtdb(rtdb,rtdb_fname,grtdbname)
      util_sgstart=.true.
      sgactive=util_sgstart
      if (.not. rtdb_open(grtdbname, 'old', rtdb)) call errquit
     *     (' rtdb_open old failed ',0,0)

      if(myproc.eq.0.and.sgprint) then
        write(6,*) ga_nodeid(),' of group',mygroup,
     *     'active = ',sgactive, 'rtdbname ',grtdbname
      endif
c
c     Create movecs for each subgroup
c
      if(myproc.eq.0.and.sgprint) then
        write (6,*) 'mygroup=',mygroup,' rtdb=',rtdb,
     *            ' grtdbname=',grtdbname
      endif
      call util_sgmovecs(rtdb)
c
c
      luout_0=luout
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
      if(sgactive) then
        if(myproc.eq.0.and.sgprint) then
          write(0,*) ga_nodeid(),' deleting cloned rtdb '
        endif
        if (.not. rtdb_close(rtdb, 'delete'))call errquit
     $        (' sgend: rtdb_close and delete failed ',0,0)
      endif
      if(mygroup.ne.1) then
        close(luout, status='keep')
      endif
      luout=luout_0
      if(myproc.eq.0) then
        write(luout,*) 'end of subgroups call '
      endif
      mypgroup_old = ga_pgroup_get_default()
      call ga_pgroup_set_default(ga_pgroup_get_world())
      ignore = ga_pgroup_destroy(mypgroup_old)
      mygroup = 1
      ngroups = 1
      sgactive = .false.
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
      write(pname,333) mygroup
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

c setup

      me=ga_nodeid()
      inode = ga_cluster_nodeid()
      nnodes = ga_cluster_nprocs(inode)
      if(ga_cluster_nnodes().gt.maxgroups) call 
     *     errquit('sggosmp: increase maxgroups',0,0) 

c make each smp node into a group

      do i=0,ga_cluster_nnodes()-1
         do j=0,ga_cluster_nprocs(i)-1
            proclist(j+1)=ga_cluster_procid(i,j)
         enddo
         grouplist(i+1)=ga_pgroup_create(proclist,ga_cluster_nprocs(i))
      enddo

      call ga_pgroup_sync(ga_pgroup_get_default())

C Set the group

      mygroup=inode+1
      write(0,'(I4,A,I4)') me,' joining group',grouplist(inode+1)
      call ga_pgroup_set_default(grouplist(inode+1))

      return
      end

      subroutine util_sginit_irreg()
      implicit none
c
c     creates custom subgroups
c
#include "stdio.fh"
#include "global.fh"
#include "util_sgroup.fh"
      integer inode,nnodes
      integer i,j,proclist(maxgroups,maxcpus), myprocid,
     $ proc_counter,group_counter
c
      inode = ga_cluster_nodeid()

C This code loops over processes and divides them into groups based 
C upon cpusperg_array() and stores this information in proclist().
C At the same time, it figures out which group the local process is in
C and stores this information in mygroup.
      proc_counter=1
      group_counter=1
      do i=0,ga_cluster_nnodes()-1
         do j=0,ga_cluster_nprocs(i)-1
           proclist(group_counter,proc_counter)=ga_cluster_procid(i,j)
           myprocid=ga_cluster_procid(inode,j)
           if (proclist(group_counter,proc_counter).eq.myprocid) 
     $          mygroup=group_counter
           if (cpusperg_array(group_counter).eq.proc_counter) then
             group_counter=group_counter+1
             proc_counter=0
           end if
           proc_counter=proc_counter+1
         enddo
      enddo

      grouplist(mygroup)=
     $  ga_pgroup_create(proclist(mygroup,1),cpusperg_array(mygroup))

      call ga_sync(ga_pgroup_get_default())
      call ga_pgroup_set_default(grouplist(mygroup))

      write (6,*) 'me=',ga_nodeid(),' is the part of subg=',mygroup

      return
      end

      subroutine util_sginit()
      implicit none
c
c     Create subgroups, this has been rewritten to work for 
c     general cases ! DW
c
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
      integer nproc,myproc,nchunkq,nremainq,npsize,i,j,n,numprocs,
     $  ncountq,meafter
      integer num_proc(0:maxgroups-1),proclist(maxgroups),
     $  proc_group(0:maxgroups-1)


c  on either smp or cluster of smp, ga_nnodes is consistent value, 
c  so ga_nnodes as starting base
      nproc = ga_nnodes()
      myproc = ga_nodeid()

c split number of procs into groups, starting group index:0
      num_proc(0)=0
      nchunkq=nproc/ngroups
      nremainq=mod(nproc,ngroups)

c for irregualr distribution, this part make sure the lower groups have 
c more procs
c for example, for 5 procs on 3 groups, 1st group  has 2 procs, 2nd has 2,
c and 3rd has 1
      if(nremainq.ne.0) then
         npsize=nchunkq+1
      else
         npsize=nchunkq
      endif

c assign chunk of procs to sub groups--the  position of the chunk 
c procs in  groups
      do i=1,ngroups-1

        if(i.lt.nremainq) then
          ncountq=nchunkq+1
        else
          ncountq=nchunkq
        end if
        num_proc(i)=npsize
        npsize=npsize+ncountq

      enddo
      num_proc(ngroups)=nproc

c assign each proc to individual group , then create individual group
      do i=1,ngroups
         do j=num_proc(i-1),num_proc(i)-1
            n=j-num_proc(i-1)+1
            proclist(n)=j
         end do
         numprocs=num_proc(i)-num_proc(i-1)
         proc_group(i-1)=ga_pgroup_create(proclist, numprocs)
      enddo

c set default individual group
      do i=1,ngroups
        if(num_proc(i-1).le.myproc.and.myproc.le.num_proc(i)-1) then
          call ga_pgroup_set_default(proc_group(i-1))
          mygroup = i
        endif
      enddo

      meafter=ga_nodeid()
      
      write(*,"(2x,'proc.',I4,' now is proc.',I4,' in group',I4)")
     $myproc,meafter,mygroup
      

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
c Create movecs.mygroup
c
       call util_pname0(movecs_in,sg_movecs_out,mygroup,ngroups)
  
      if(myproc.eq.0.and.sgprint) then
        write(6,*) 'mygroup=',mygroup,
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
      util_sgroup_mygroup = mygroup
      end
      integer function util_sgroup_numgroups()
#include "util_sgroup.fh"
      util_sgroup_numgroups = ngroups
      end

