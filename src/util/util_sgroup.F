C$Id: util_sgroup.F,v 1.10 2007-02-14 17:10:49 windus Exp $
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
      logical lstatus
c
c Read number of subgroups from input
c
      if (.not.rtdb_get(rtdb, 'subgroups_number', mt_int, 1,
     &   ngroups))ngroups=1
      gosmp=.false.
      if (ngroups.eq.ga_cluster_nnodes()) then
       gosmp=.true.
       ngroups=ga_cluster_nnodes()
      endif
      if (ngroups.gt.ga_cluster_nnodes()) then
       call errquit('Too many subgroups are requested',0,0)
      endif

      write (6,*) 'ngroups=',ngroups,' gosmp=',gosmp
c
c     preliminary rtdb cloning
c
      if (.not. rtdb_getfname(rtdb, rtdb_fname)) call 
     *     errquit('rtdb_fname call failed',0,0)
      if (.not. rtdb_close(rtdb, 'keep')) call errquit(
     *     ' failed to close rtdb for sg ',0,0)
c
c    Create groups
c   
      write(luout,*) 'Creating groups'
      if(gosmp) then
         call util_sginit()
         call util_sgrtdb(rtdb,rtdb_fname,grtdbname)
         util_sgstart=.true.
      else
c        call util_sginit_irreg()
         call util_sginit()
         call util_sgrtdb(rtdb,rtdb_fname,grtdbname)
         util_sgstart=.true.
      endif
      sgactive=util_sgstart
      if (.not. rtdb_open(grtdbname, 'old', rtdb)) call errquit
     (     (' rtdb_open old failed ',0,0)
      write(6,*) ga_nodeid(),' of group',mygroup,
     *     'active = ',sgactive, 'rtdbname ',grtdbname
c
c     Create movecs for each subgroup
c
      write (6,*) 'mygroup=',mygroup,' rtdb=',rtdb,
     $            ' grtdbname=',grtdbname
      call util_sgmovecs(rtdb)
c
c     nasty hack stdoutput redirect to fort.90+group for 
c     groupid gt 0
c
      luout_0=luout
c      if(mygroup.ne.1.and.gosmp) then
c        luout=90+mygroup
c        write(6,*) mygroup,'new value for luout',luout
c        if(ga_nodeid().eq.0) then
c           inquire(luout, exist=lstatus,iostat=errout)
c           write(6,*) 'lstatus=',lstatus
c           if(lstatus) then
c              open(luout,  status='old',iostat=errout)
c              if(errout.eq.0)close(luout, status='delete')
c           else
c              write(6,*) ' fort.9x not there'
c           endif
c           open(luout, status='new')
c           write(6,*) 'after the open statement'
c        endif
c        call ga_sync()
c      endif
      write(6,*) 'everyone should have an open file at this point'
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
      if(mygroup.ne.1) then
        close(luout, status='keep')
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
      character*256 scratchrtdb,cmd,cmd1
      character*10 pname
      integer util_system
      external util_system
c
      write(luout,*) 'util_sgrtdb start'
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
               Write (6,*) 'cmd done'
            else
               Write (0,*) cmd, ' failed '
               call errquit('util_sgrtdb: problem with system call',
     *              0, 0)
            endif
         endif
         grtdbname=scratchrtdb
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

      subroutine util_sginit_irreg()
      implicit none
c
c     creates custom subgroups
c
#include "stdio.fh"
#include "global.fh"
#include "util_sgroup.fh"
      integer me,inode,nnodes
      integer i,j,proclist(maxgroups,maxcpus), myprocid,
     $ proc_counter,group_counter
c
      me=ga_nodeid()
      inode = ga_cluster_nodeid()
!      nnodes = ga_cluster_nprocs(inode)

      write(6,*) 'me in util_sginit_irreg is ', me
      proc_counter=1
      group_counter=1
      write(6,*) 'ga_cluster_nnodes ',ga_cluster_nnodes()
      do i=0,ga_cluster_nnodes()-1
         write(6,*) 'ga_cluster_nprocs(i) is ',ga_cluster_nprocs(i),i
         do j=0,ga_cluster_nprocs(i)-1
           write(6,*) 'i and j ', i, j
           write(6,*) 'group_counter,proc_counter,i,j,ga_clus_procid'
     $       , group_counter,proc_counter,i,j,ga_cluster_procid(i,j)
           proclist(group_counter,proc_counter)=
     $     ga_cluster_procid(i,j)
           myprocid=ga_cluster_procid(inode,j)
           write(6,*) 'myprocid,proclist ',myprocid,
     $        proclist(group_counter,proc_counter),
     &        group_counter,proc_counter
c          if (proclist(i,j).eq.myprocid) mygroup=group_counter
           if (proclist(group_counter,proc_counter).eq.myprocid) 
     $          mygroup=group_counter

           if (cpusperg_array(group_counter).eq.proc_counter) then
             group_counter=group_counter+1
             proc_counter=0
           end if
           proc_counter=proc_counter+1
         enddo
      enddo

      write(6,*) 'mygroup,cpusperg_array ',mygroup,
     $     cpusperg_array(mygroup)
      grouplist(mygroup)=
     $ga_pgroup_create(proclist(mygroup,1),cpusperg_array(mygroup))

      call ga_sync()
      call ga_pgroup_set_default(grouplist(mygroup))

      me=ga_nodeid()
      write (6,*) 'me=',me,' is the part of subg=',mygroup

      return
      end

      subroutine util_sginit()
      implicit none
c
c     Create subgroups
c
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
      integer cur_mygroup,def_group,me

      def_group = ga_pgroup_get_default()
      cur_mygroup=ga_pgroup_split(def_group,ngroups)
      if (cur_mygroup.eq.-1) then
         call errquit('Died in ga_pgroup_split',0,0)
      endif
      call ga_pgroup_set_default(cur_mygroup)
      mygroup=ga_pgroup_get_default()

      me=ga_nodeid()
      write (6,*) 'me=',me,' is the part of subg=',mygroup

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
c
c Tested only for the case when movecs are created
c

c
c Read movecs filename with path
c
       if (.not.(rtdb_cget(rtdb,'scf:input vectors',1,movecs_in).or.
     $           rtdb_cget(rtdb,'dft:input vectors',1, movecs_in)))
     $ movecs_in = 'atomic'
        if (movecs_in.eq.'atomic')
     $  call util_file_name('movecs',.true.,.false.,movecs_in)
c
c Create movecs.mygroup
c
        call util_pname0(movecs_in,sg_movecs_out,mygroup,ngroups)
  
        write(6,*) 'mygroup=',mygroup,
     $             ' movecs_in',movecs_in,
     $             ' sg_movecs_out=',sg_movecs_out
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
      if (.not. rtdb_cput(rtdb,tag,1,
     $    sg_movecs_out)) then
       write (6,*) 'error in rtdb_cput'
      endif
c
      return
      end

