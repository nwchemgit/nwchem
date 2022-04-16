C$Id$
C>
C> \defgroup subgrp NWChem interface to GA sub-groups
C>
C> \ingroup subgrp
C> @{
C>
C> \file util_sgroup.F
C> \brief The NWChem interface to the GA sub-groups 
C>
C> This file contains routines that facilitate the use of GA sub-groups
C> in NWChem [1-3]. The routines take care of issues such as,
C> redirecting the standard output from different sub-groups,
C> replicating the runtime database, accessing the movecs file, etc. In
C> addition it provides mechanisms for creating pools of processors in
C> a number of different ways. 
C>
C> To date there is no input block to control this functionality.
C> Instead set directives that put data directly into the runtime
C> database are required.
C>
C> ### References ###
C>
C> [1] T. L. Windus, S. M. Kathmann, L. D. Crosby,
C>     "High performance computations using dynamical nucleation
C>     theory", <i>J. Phys.: Conf. Ser.</i> (2008) <b>125</b>,
C>     012017, DOI:
C>     <a href="https://doi.org/10.1088/1742-6596/125/1/012017">
C>     10.1088/1742-6596/125/1/012017</a>.
C>
C> [2] B. J. Palmer, S. M. Kathmann, M. Krishnan, V. Tipparaju,
C>     J. Nieplocha, "The use of processor groups in molecular
C>     dynamics simulations to sample free-energy states",
C>     <i>J. Chem. Theory Comput.</i> (2007) <b>3</b>, 583-592,
C>     DOI:
C>     <a href="https://doi.org/10.1021/ct600260u">
C>     10.1021/ct600260u</a>.
C>
C> [3] J. Nieplocha, B. Palmer, V. Tipparaju, M. Krishnan, H. Trease,
C>     E. Apra, "Advances, applications and performance of the Global
C>     Arrays shared memory programming toolkit", 
C>     <i>Int. J. High Perf. Comput. Appl.</i> (2006) <b>20</b>,
C>     203-231, DOI:
C>     <a href="https://doi.org/10.1177/1094342006064503">
C>     10.1177/1094342006064503</a>.
C>
      block data groups
#include "util_sgroup.fh"
C     Make sure that the group code works even if groups are not set up yet
      data mygroup(0) /1/
      data my_ga_grp(0)/-1/
      data my_ga_0_grp(0)/-999/ ! Bogus
      data ngroups(0) /1/
      data sgprint /.true./
      data depth /0/
      end

C###############################################################################
C>
C> \brief Initialize the subgroup layer
C>
C> \return Returns .TRUE. if successful (does not return otherwise)
C>
      logical function util_sgstart(rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb !< [Input] The RTDB handle
      integer groups_want
      integer dum(1)
      integer method
C added option to allow writting of group rtdb's
C to permanent_dir(.false.) or scratch_dir(.true.)
      integer dir
c
      dir = 1 ! permanent_dir by default
c
c Read number of subgroups from input
c
      if (.not.rtdb_get(rtdb, 'subgroups_number', mt_int, 1,
     &   groups_want))groups_want=1

      method = 1
      call util_sggo(rtdb,groups_want,method,dum,dir)

      util_sgstart=.true.
      return
      end
C###############################################################################
C>
C> \brief Do the actual initialization
C>
      subroutine util_sggo(rtdb,groups_want,method,array_cpu,dir)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "errquit.fh"

c
c Input
c
      integer rtdb         !< [Input] The RTDB handle
      integer groups_want  !< [Input] The number of subgroups requested
      integer method       !< [Input] The method to partition the 
                           !< processor pool
                           !< * 1 -- use `groups_want` to generate equal
                           !<   sized groups (`array_cpu` ignored)
                           !< * 2 -- turn each SMP box into a group
                           !<   (`array_cpu` and `groups_want` ignored)
                           !< * 3 -- use `array_cpu(groups_want)` to
                           !<   define number of nodes per group
                           !< * 4 -- use `array_cpu(groups_want+nnodes)`
                           !<   to define which nodes per group
      integer array_cpu(*) !< [Input] Specification of the subgroups
c
c     rtdb is the runtime database file
c
c     method is:
c     1 -- use groups_want to generate equal sized groups (array_cpu ignored)
c     2 -- turn each SMP box into a group (array_cpu and groups_want ignored)
c     3 -- use array_cpu(groups_want) to define number of nodes per group
c     4 -- use array_cpu(groups_want+nnodes) to define which nodes per group
c
      integer i
      integer meinworld, meafter
      integer X0(1)
      integer idum

C strings for fname definition
      character*256 permdir
      character*256 fprefix
      integer myproc,nproc
      integer mypgroup_old
C
C     Added option to allow group rtdb's to be written
C     to permanent_dir (0) or scratch_dir(1)
      integer dir
      logical dir2
#ifndef GANXTVAL
C  Turn off GA load balancing in groups
      logical util_statldb
      external util_statldb
      logical dumlog
#endif

C     convert to logical
      if (dir .gt. 0) then
        dir2 = .true.
      else
        dir2 = .false.
      endif
#ifndef GANXTVAL
C     Turn off dynamic load balancing, since it fails in groups
C     We never turn this back on
      if (groups_want .ne. 1) dumlog = util_statldb(.true.,rtdb)
#endif
c
c     dont use disk for gridpts, etc. in DFT
c     We never turn this back on
c
      if(.not.rtdb_put(rtdb,'dft:largenode',mt_log,1,.true.))
     &   call errquit('util_sggo: rtdb_put nodisk failed',1,RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:noio', mt_int, 1, 1))
     &   call errquit('util_sggo: rtdb_put noio failed',1,RTDB_ERR)
      if(.not.rtdb_put(rtdb,'dft:xcreplicated',mt_log,1,.false.))
     &  call errquit('util_sggo: rtdb_put xcreplicated fail',1,RTDB_ERR)

c

      mypgroup_old = ga_pgroup_get_default()
      if (mypgroup_old .eq. ga_pgroup_get_world()) then
        depth = 0
        mygroup(depth) = 1
        ngroups(depth) = 1
        my_ga_grp(depth) = ga_pgroup_get_world()
        x0(1) = 0
        if (my_ga_0_grp(0) .eq. -999) then
          my_ga_0_grp(depth) = ga_pgroup_create(x0,1)
        end if

c       preliminary rtdb cloning
C       Work around the 36 char cutoff for fname
        if (.not. rtdb_getfname(rtdb, rtdb_fname(depth))) call
     *     errquit('rtdb_fname call failed',0,0)
        write(permdir,'(256(a))') (' ', i=1,256)
        write(fprefix,'(256(a))') (' ', i=1,256)
        if(.not. rtdb_cget(rtdb,'permanent_dir',1,permdir)) then
           permdir = '.'
        endif
        if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix))
     &     call errquit('rtdb get file_prefix failed',0,0)
        idum = 1
  747   if(permdir(idum:idum).eq. ' ' .and. idum.lt.inp_strlen(permdir))
     &                then
          idum = idum + 1
          goto 747
        endif
        write(rtdb_fname(depth),'(a,a,a,a)')
     &    permdir(idum:inp_strlen(permdir)),
     &    '/',
     &    fprefix(1:inp_strlen(fprefix)),
     &    '.db'
      end if
      if (.not. rtdb_close(rtdb, 'keep')) call errquit(
     *     ' failed to close rtdb for sg ',0,0)

      depth = depth + 1
      if (depth .gt. maxdeep) 
     *     call errquit('depth is greater than maxdeep',0,0)

      myproc=ga_nodeid()
      nproc=ga_nnodes()
c      sgprint = util_print('sgroup info', print_debug)
c
c    Create groups
c   
      if(myproc.eq.0) write(LuOut,*) 'Creating groups'

c initialize sub groups using a method

      if (method .eq. 1) then
        call util_sginit(groups_want)
      else if (method .eq. 2) then
        call util_sginit_smp()
      else if (method .eq. 3) then
        call util_sginit_irreg(array_cpu,groups_want)
      else if (method .eq. 4) then
        call util_sginit_power(array_cpu,groups_want+nproc)
      else
        call errquit('util_sggo: bad method',0,0)
      endif

      if (sgprint) then
        meafter=ga_nodeid()
        meinworld = ga_pgroup_nodeid(ga_pgroup_get_world())
        write(LuOut,
     *   "(' proc ',I4,' was (',
     *       I4,' ',I4,' ',I4,
     *       ') now is (',
     *       I4,' ',I4,' ',I4,
     *       ')')")
     *   meinworld,
     *   myproc,mygroup(depth-1),my_ga_grp(depth-1),
     *   meafter,mygroup(depth),my_ga_grp(depth)
         call util_flush(LuOut)
      endif

      call ga_pgroup_sync(mypgroup_old) ! Make sure all is well

      call util_sgrtdb(rtdb,rtdb_fname(depth-1),rtdb_fname(depth),dir2)
      if (.not. rtdb_open(rtdb_fname(depth), 'old', rtdb)) call errquit
     *     (' rtdb_open old failed ',0,0)

      if(myproc.eq.0.and.sgprint) then
       write(LuOut,*) ga_nodeid(),' of group',mygroup(depth),'rtdbname '
     *           ,rtdb_fname(depth), ' rtdb=',rtdb
      endif
c
c     Create movecs for each subgroup
c
      call util_sgmovecs(rtdb)
c
c
      if(myproc.eq.0.and.sgprint) then
        write(LuOut,*) 'everyone should have an open file at this point'
        call util_flush(LuOut)
      endif
c     These next lines are since we are now done creating groups
      call util_flush(LuErr)
      call util_flush(LuOut)
      call ga_pgroup_sync(mypgroup_old) ! Make sure all is well
      return
      end
C###############################################################################
C>
C> \brief Finalize the subgroup layer
C>
      subroutine util_sgend(rtdb)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"
      integer rtdb   !< [Input] The RTDB handle
      integer myproc
      logical ignore
      integer mypgroup_old

      mypgroup_old = ga_pgroup_get_default()
      myproc=ga_nodeid()
      if (mypgroup_old .eq. ga_pgroup_get_world()) then
        write(LuErr,*) ' Node ',myproc,' tried to close world group '
        return
      endif
      myproc=ga_nodeid()
      if(myproc.eq.0.and.sgprint) then
        write(LuErr,*) ' deleting cloned rtdb '
      endif
C     Make sure all nodes are done using RTDB file
      call ga_pgroup_sync(mypgroup_old)
      if (.not. rtdb_close(rtdb, 'delete'))call errquit
     $        (' sgend: rtdb_close and delete failed ',0,0)
      if(myproc.eq.0) then
        write(LuOut,*) 'Closing subgroup '
      endif
      if (sgprint) call util_flush(LuOut)

      call ga_pgroup_set_default(my_ga_grp(depth-1))
      if (myproc .eq. 0) then 
      ignore = ga_pgroup_destroy(my_ga_0_grp(depth))
      endif
      ignore = ga_pgroup_destroy(mypgroup_old)
      my_ga_0_grp(depth) = -999
      my_ga_grp(depth) = -999
c
c     reopen higher up rtdb
c
      depth = depth - 1
      if (.not. rtdb_open(rtdb_fname(depth), 'old', rtdb))
     $     call errquit('sgend: rtdb_open old failed ', 0, 0)

      return
      end

C###############################################################################
C>
C> \brief Copy the RTDB for each subgroup
C>
C> Copy the RTDB for each subgroup. This way the subgroups can 
C> independently modify the RTDB without messing each other's work up.
C> Depending on the kind of calculation one may want to store the 
C> RTDB replicas in the scratch directory or the permanent directory.
C> This behavior can be controlled by the `dir` argument.
C>
      subroutine util_sgrtdb(rtdb,name_in,name_out,dir)
      implicit none
      integer rtdb           !< [Input] The RTDB handle
      character*256 name_in  !< [Input] The filename of the global RTDB
      character*256 name_out !< [Output] The filename of the subgroup
                             !< RTDB
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util_sgroup.fh"
      integer ii, iend, ibegin
      logical ortdbstate,ortdbstate2,status
      character*256 grouprtdb
      character*256 groupname
      integer myproc

C     Added option to allow group rtdb's to be written to
C     Permanent_dir or scratch_dir
      logical dir !< [Input] Flag determining the location of the 
                  !< replica RTDBs
                  !< * .FALSE. changes to `permanent_dir`
                  !< * .TRUE. changes to `scratch_dir`

      myproc=ga_nodeid()

      if(myproc.eq.0.and.sgprint) then
        write(LuOut,*) 'util_sgrtdb start'
        call util_flush(LuOut)
      endif
c
c     do a copy of rtdb to permanent_dir or scratch_dir for each group
c
      ibegin = 1
      iend   = ibegin + (6-1)
      write(groupname(ibegin:iend),'(6A)') 'sg_db.'
      do ii=1,depth
        ibegin = iend + 1
        iend   = ibegin + (4-1)
C       9999 is largest 4 digit number, and 256 is just plain long (used above)
C       Of course a long group name will most likely fail long before
C       it gets to that long - other limitations on file name length exist
        if (iend .gt. 256 .or. mygroup(ii) .gt. 9999)
     *    call errquit(' util_sgrtdb: groups too deep',0,0)
        write(groupname(ibegin:iend),'(i4.4)') mygroup(ii)
      end do

      call util_file_name(groupname(1:iend),dir,.false.,grouprtdb)
C
C     Do not want to have other nodes to be using old RTDB file anymore
      call ga_pgroup_sync(ga_pgroup_get_default())
C
C     Node 0 of the group now copies the original RTDB to a group specific one
C
      if (myproc .eq. 0) then
        if(sgprint) then
           write(LuOut,*)'world rtdbname is ',
     &       name_in(1:inp_strlen(name_in))
           write(LuOut,*)'group rtdbname is ',
     &       grouprtdb(1:inp_strlen(grouprtdb))
           call util_flush(LuOut)
        endif
        inquire(file=name_in(1:inp_strlen(name_in)),exist=status)
        if(.not.status) then
c          most likely cause of error: permdir not available on this node
           write(LuErr,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(LuErr,*) ' ! please change permanent_dir to a   !'
           write(LuErr,*) ' ! filesystem available to all the    !'
           write(LuErr,*) ' ! processors (e.g. NFS, GPFS, Lustre)!'
           write(LuErr,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           call errquit(' util_sgrtdb: failed to open rtdb',0,0)
        endif
        call util_file_copy(name_in(1:inp_strlen(name_in)),
     &        grouprtdb(1:inp_strlen(grouprtdb)))
        inquire(file=grouprtdb(1:inp_strlen(grouprtdb)),
     &          exist=status)
        if (.not. status) then
          write (LuErr,*) 'copy failed'
          call errquit('util_sgrtdb: problem with system call', 0, 0)
        else if (sgprint) then
          write(LuOut,*) ' rtdb copied ', ga_pgroup_get_default()
          call util_flush(LuOut)
        endif
      endif
      name_out=grouprtdb
C     Do not want to have other nodes using new RTDB before it exists
      call ga_pgroup_sync(ga_pgroup_get_default())
      return
      end
C###############################################################################
C>
C> \brief Initialize the subgroup layer using the SMP method
C>
C> The SMP method targets machines that are clusters of shared memory
C> nodes. It creates subgroups where each subgroup coincides with
C> a shared memory node.
C>
      subroutine util_sginit_smp()
      implicit none
c
c     creates one subgroup for each smp
c
#include "global.fh"
#include "util_sgroup.fh"
      integer i,j,proclist(maxcpus),num_proc

c setup
      if (depth .ne. 1) then
        call errquit(' util_sginit_smp: bad depth',0,0)
      endif
      ngroups(depth) = ga_cluster_nnodes()
      mygroup(depth) = 1 + ga_cluster_nodeid()
      my_ga_grp(depth) = -99

c make each smp node into a group

      do i=0,ngroups(depth)-1
        num_proc = ga_cluster_nprocs(i)
        if (num_proc .gt. maxcpus) call errquit('increase maxcpus',0,0)
        do j=0,num_proc - 1
          proclist(j+1)=ga_cluster_procid(i,j)
        enddo
        if (i+1 .eq. mygroup(depth)) then
          my_ga_grp(depth) = ga_pgroup_create(proclist,num_proc)
        endif
      enddo

c create the zero group
     
      num_proc = ngroups(depth)
      if(num_proc .gt. maxcpus)call errquit('increase maxcpus',0,0)
      do i=0, num_proc - 1
        proclist(i+1)=ga_cluster_procid(i,0)
      enddo
      my_ga_0_grp(depth) = ga_pgroup_create(proclist,num_proc)

C Set the group

      call ga_pgroup_set_default(my_ga_grp(depth))

      return
      end
C###############################################################################
C>
C> \brief Initialize the subgroup layer with irregularly sized groups
C>
C> This method generates irregularly sized subgroups in a simple way.
C> The array `cpusperg_array` (CPUs per Group) specifies the size of
C> every subgroup. This routine generates subgroups to match that 
C> specification.
C>
      subroutine util_sginit_irreg(cpusperg_array,length)
      implicit none
c
c     creates custom subgroups using cpusperg_array()
c
#include "global.fh"
#include "util_sgroup.fh"
      integer nproc, myproc
      integer i,j,proclist(maxcpus), proc_counter,group_counter
      integer length                 !< [Input] The number of subgroups
      integer cpusperg_array(length) !< [Input] The size of each
                                     !< subgroup

c setup

      nproc = ga_nnodes()
      myproc = ga_nodeid()
      my_ga_grp(depth) = -99
      mygroup(depth)=-99
      ngroups(depth)=length

C verify input all at once instead of constantly in loops

      proc_counter = 0
      do i=1,length
        if(cpusperg_array(i) .le. 0) then
          call errquit(' util_sginit_irreg: zero cpusperg_array',0,0)
        endif
        if(cpusperg_array(i) .gt. maxcpus) then
          call errquit(' util_sginit_irreg:increase maxcpus',0,0)
        endif
        proc_counter = proc_counter + cpusperg_array(i)
      enddo
      if (proc_counter .ne. nproc) then
        call errquit(' util_sginit_irreg:bad cpusperg_array',0,0)
      end if

C make groups
      proc_counter=1
      group_counter=1
      do i=0,nproc-1
        proclist(proc_counter) = i
        if (i.eq.myproc) then
          mygroup(depth)=group_counter
        endif
        if (cpusperg_array(group_counter).eq.proc_counter) then
          if(mygroup(depth) .eq. group_counter) then
            my_ga_grp(depth)=ga_pgroup_create(proclist,proc_counter)
          endif 
          group_counter=group_counter+1
          proc_counter=0
        end if
        proc_counter=proc_counter+1
      enddo

c create the zero group

      if(ngroups(depth).gt.maxcpus)call errquit('increase maxcpus',0,0)
      proc_counter=1
      group_counter=1
      do i=0,nproc-1
        if (proc_counter .eq. 1) then
          proclist(group_counter) = i
        end if
        if (cpusperg_array(group_counter).eq.proc_counter) then
          group_counter=group_counter+1
          proc_counter=0
        end if
        proc_counter=proc_counter+1
      enddo
      my_ga_0_grp(depth) = ga_pgroup_create(proclist,ngroups(depth))

C Set the group

      call ga_pgroup_set_default(my_ga_grp(depth))

      return
      end
C###############################################################################
C>
C> \brief Initialize the subgroup layer using the "power" method
C>
C> The power method creates subgroups based on the specification in
C> `cpu_array`. The data in this array specifies the size of each
C> subgroup as well as which processes are part of each subgroup.
C> The data is stored in a vector of the following form
C> \f{eqnarray*}{
C>   \left(G_1,n_1,n_2,n_3,\ldots,n_{G_1},G_2,n_1,n_2,n_3,\ldots,
C>   n_{G_2},\ldots\right)
C> \f}
C> where \f$G_1\f$, \f$G_2\f$, etc. are the sizes of the respective
C> subgroups, and \f$n_1, \ldots, n_G\f$ are the processor ranks that
C> are part of that subgroup.
C>
      subroutine util_sginit_power(cpu_array,length)
      implicit none
c
c     creates custom highly custom subgroups using cpu_array()
c     is is of the form G1,n1,n2,n3,...,nn,G2,n1,n2,n3,...,nn,...
c     Where the G's are the size of the groups, and the n's
c     are the specific nodes in the groups
c     length = nnodes + ngroups
c
#include "global.fh"
#include "util_sgroup.fh"
      integer nproc, myproc
      integer i,j,proclist(maxcpus), proc_counter,group_counter
      integer pos_counter, group_size
      integer length, cpu_array(length)

c setup

      nproc = ga_nnodes()
      myproc = ga_nodeid()
      my_ga_grp(depth) = -99
      mygroup(depth)=-99
      ngroups(depth)=-99

C make groups
      group_counter=0
      pos_counter=0
      do while (pos_counter .lt. length)
        group_counter = group_counter + 1
        pos_counter = pos_counter + 1
        group_size = cpu_array(pos_counter)
        proc_counter = 0
        do i=1,group_size
           pos_counter = pos_counter + 1
           proc_counter=proc_counter+1
           if (proc_counter .gt. maxcpus) then
             call errquit('increase maxcpus',0,0)
           end if
           proclist(proc_counter) = cpu_array(pos_counter)
           if (proclist(proc_counter) .eq. myproc) then
             mygroup(depth)=group_counter
           endif
        end do
        if(mygroup(depth) .eq. group_counter) then
            my_ga_grp(depth)=ga_pgroup_create(proclist,proc_counter)
        endif
      end do
      ngroups(depth)=group_counter

c create the zero group

      if(ngroups(depth).gt.maxcpus)call errquit('increase maxcpus',0,0)
      group_counter=0
      pos_counter=1
      do while (pos_counter .lt. length)
        group_counter = group_counter + 1
        group_size = cpu_array(pos_counter)
        pos_counter = pos_counter + 1
        proclist(group_counter) = cpu_array(pos_counter)
        pos_counter = pos_counter + group_size
      end do
      my_ga_0_grp(depth) = ga_pgroup_create(proclist,ngroups(depth))

C Set the group

      call ga_pgroup_set_default(my_ga_grp(depth))

      return
      end


C###############################################################################
C>
C> \brief Initialize the subgroup layer with regular subgroup sizes
C>
      subroutine util_sginit(groups_want)
      implicit none
#include "global.fh"
#include "util_sgroup.fh"
c
c     Create subgroups of a constant size
c
      integer groups_want !< [Input] The number of subgroups wanted
      integer nproc,myproc,nchunkq,nremainq,i,j,n
      integer num_proc,proclist(maxcpus)

C setup

      ngroups(depth) = groups_want
      nproc = ga_nnodes()
      myproc = ga_nodeid()
      mygroup(depth) = -99
      my_ga_grp(depth) = -99

      if (nproc .lt. ngroups(depth)) then 
         if(myproc .eq. 0) write(*,*) ' Will only create group of',nproc
         ngroups(depth) = nproc
      endif
      nchunkq=nproc/ngroups(depth)
      nremainq=mod(nproc,ngroups(depth))
      if ((nchunkq .gt. maxcpus) .or.
     *     (nchunkq+1 .gt. maxcpus .and. nremainq .ne. 0)) then
        call errquit('increase maxcpus',0,0)
      end if

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
         if (mygroup(depth) .eq. i) then
            my_ga_grp(depth) = ga_pgroup_create(proclist, num_proc)
         endif
      enddo

c create the zero group

      if(ngroups(depth).gt.maxcpus)call errquit('increase maxcpus',0,0)
      n = 0
      do i=1,ngroups(depth)
         num_proc = nchunkq
         if(i .le. nremainq) num_proc = num_proc + 1
         proclist(i) = n
         n = n + num_proc
      enddo
      my_ga_0_grp(depth) = ga_pgroup_create(proclist,ngroups(depth))

c set the group

      call ga_pgroup_set_default(my_ga_grp(depth))

      return
      end
C###############################################################################
C>
C> \brief Create movecs file names for the current subgroup
C>
C> The movecs file holds the SCF vectors of a particular calculation
C> and is therefore a crucial data object. With subgroups every group
C> works on a calculation of its own and therefore needs a movecs file
C> of its own. This routine creates and stores the name of the movecs
C> file for the current subgroup.
C>
      subroutine util_sgmovecs(rtdb)
      implicit none
      integer rtdb !< [Input] The RTDB handle
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
C
c Read movecs filename with path
c
      myproc=ga_nodeid()

       if (.not.(rtdb_cget(rtdb,'scf:input vectors',1,movecs_in).or.
     >           rtdb_cget(rtdb,'pspw:input vectors',1,movecs_in).or.
     $           rtdb_cget(rtdb,'dft:input vectors',1, movecs_in)))
     $     movecs_in = 'atomic'
       if (movecs_in.eq.'atomic')
     $     call util_file_name('movecs',.true.,.false.,movecs_in)
c
c Create movecs.mygroup(depth)
C This just keeps appending to the end
c
      call util_pname0(movecs_in,sg_movecs_out,mygroup(depth),
     *     ngroups(depth))
  
      if(myproc.eq.0.and.sgprint) then
        write(LuOut,*) 'mygroup=',mygroup(depth),
     $             ' movecs_in',movecs_in,
     $             ' sg_movecs_out=',sg_movecs_out
      endif
c
c Set output movecs
c
c      if (.not.(rtdb_cget(rtdb,'task:theory',1,theory)))
c     *     call errquit('rtdb_get task:theory failed',0,0)
c      if (theory.eq.'dft') then
c         tag='dft:output vectors'
c      else
c         tag='scf:output vectors'
c      endif
c
      tag='dft:output vectors'
      if (.not. rtdb_cput(rtdb,tag,1,sg_movecs_out)) then
       write (LuOut,*) 'error in rtdb_cput'
      endif
      tag='scf:output vectors'
      if (.not. rtdb_cput(rtdb,tag,1,sg_movecs_out)) then
       write (LuOut,*) 'error in rtdb_cput'
      endif

      tag='pspw:input vectors'
      if (.not. rtdb_cput(rtdb,tag,1,sg_movecs_out)) then
       write (LuOut,*) 'error in rtdb_cput'
      endif
      tag='pspw:output vectors'
      if (.not. rtdb_cput(rtdb,tag,1,sg_movecs_out)) then
       write (LuOut,*) 'error in rtdb_cput'
      endif

c
      return
      end
c     Hokey functions to keep others away from my common blocks
C###############################################################################
C>
C> \brief Obtain my subgroup number
C>
      integer function util_sgroup_mygroup()
#include "util_sgroup.fh"
      util_sgroup_mygroup = mygroup(depth)
      return
      end
C###############################################################################
C>
C> \brief Obtain the total number of subgroups
C>
      integer function util_sgroup_numgroups()
#include "util_sgroup.fh"
      util_sgroup_numgroups = ngroups(depth)
      return
      end
C###############################################################################
C>
C> \brief Obtain the subgroup of all rank 0 processes in each subgroup
C>
      integer function util_sgroup_zero_group()
#include "util_sgroup.fh"
      util_sgroup_zero_group = my_ga_0_grp(depth)
      return
      end
C###############################################################################
C>
C> \brief Obtain the current subgroup nesting depth
C>
      integer function util_sgroup_depth()
#include "util_sgroup.fh"
      util_sgroup_depth = depth
      return
      end
C###############################################################################
C> 
C> \brief Set the file name for standard output for this subgroup
C>
      integer function util_sgroup_set_ioname(filename)
      implicit none
      character*(*) filename !< [Input] The base file name

#include "stdio.fh"
#include "util_sgroup.fh"

      logical       found
      character*255 full_filename
      character ch_tmp

#ifdef USE_SUBGROUPS
      LuOut = mygroup(depth)+200

      call util_file_name_noprefix(filename,.false.,.false.,
     >                             full_filename)

      inquire(file=full_filename,exist=found)
*     **** FILE already exists - parse to EOF ****
      if (found) then
         open(unit=LuOut,file=full_filename,form='formatted',
     >             status='old')
         do while(.true.)
            read(LuOut,*,ERR=30,END=30) ch_tmp
         end do
 30      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1)
         backspace LuOut
#endif
*     **** FILE does not exist ****
      else
         open(unit=LuOut,file=full_filename,form='formatted')
      end if
      util_sgroup_set_ioname = LuOut
#else
      util_sgroup_set_ioname = -1
#endif
      return
      end
C###############################################################################
C>
C> \brief Switch standard output back to the default Fortran scenario
C>
      integer function util_sgroup_unset_io()
      implicit none
#include "stdio.fh"
#include "util_sgroup.fh"
      
#ifdef USE_SUBGROUPS
      close(LuOut)
      LuOut    = 6
      util_sgroup_unset_io = LuOut
#else
      util_sgroup_unset_io = -1
#endif
      return
      end
C###############################################################################
C>
C> @}
