      function mc_main(rtdb,grad,thr)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "util_sgroup.fh"  !for subgroup support
#include "const_data.fh"
#include "subgr.fh"
      logical mc_main
      logical md_driver,mc_driver,mc_init
      integer mcsteps
      logical grad
      external grad

      logical status
      integer natom 
      integer c_start
      integer vel_init,acc_init
      logical forward,backward   ! tells you which part of the traj is being called
      logical mc_data_set_ifirc, mc_data_get_ifirc
      integer mc_data_get_forside, mc_data_get_forside0
      integer mc_data_get_backside, mc_data_get_backside0
      integer i,j,k,naccept 
      integer c_array, vel, acc, prp 
      integer mc_data_get_i_c_array, md_data_get_i_v, md_data_get_i_a 
      integer md_data_get_i_c, mc_data_get_i_in_vel,mc_data_get_i_in_acc
      integer mc_data_get_natom, mc_data_get_i_prp
      integer tag,charge,md_data_get_i_t,md_data_get_i_q
      character*32 thr
      logical ircflag,mc_data_set_forward,mc_data_set_backward
      integer i_sc, mc_data_get_i_sc, nonreact,reactive,reverse
      integer i_s, mc_data_get_i_s,mc_data_get_mcsteps
      double precision avg,rc,temp,mc_data_get_temp
      logical md_data_set_temp
      logical mc_data_set_side_prev, mc_data_set_nxing
      logical mc_data_set_trajnum
      double precision mdtemp,mdtemp10, md_data_get_temp
      logical dbug

      integer g_a,ncpu
C     PROCESS AND SUBGROUP VARS
      integer myid,idbig,inodesbig
      integer idmedium, idzero, g_pr, numgroups
      integer l_groupnums, k_groupnums
C   I/O vars
      character*4 prcfil
      character*4 prcfil2
      character*256 fprefix, fprefix2
      character*256 mc_dir
      integer iw 
      iw=0
      if (ga_nodeid().eq.0) then
       iw=6
      endif
      call ga_sync()
C     INFO ON PROCESSES in the big group
      myid = ga_nodeid()
      idbig = ga_pgroup_get_default()
      inodesbig = ga_nnodes() 
c      return
      if (.not. rtdb_get(rtdb,'subgroups_number', MT_INT, 1, numgroups))
     & numgroups = ga_nnodes() 
      write(iw,*) "numgroups is ",numgroups
C in general I want nproc/want_groups
C at the moment I want numgroups=nproc
C     Allocate memory for array of group numbers
      if (.not. ma_push_get(MT_INT, numgroups, 'mc: groupnums',
     &    l_groupnums, k_groupnums)) 
     &    call errquit('mc_allocate_arrays: error groupnums',0,MA_ERR)
      write(iw,*) 'SPLITTING INTO SUBGROUPS'
C    SPLITTING INTO SUBGROUPS

C   Close the output file and reopen it as series of ".out" files
c      flush(6)  
      call ga_sync()
      CLOSE(UNIT = 6)            
      call mc_setup_ga
     &    (rtdb, myid, idbig, inodesbig,
     &    idmedium, idzero, g_pr, int_mb(k_groupnums),numgroups) ! Now in Subgroups

      write(iw,*) 'DONE SPLITTING INTO SUBGROUPS'
C Initialize group related information
C mostly file opening for each group,etc. 
      if (ga_nodeid() .eq. 0) then
         call mc_setup_group  !deals with file i/o
     &     (rtdb, idmedium, myid, int_mb(k_groupnums),prcfil, prcfil2,
     &     fprefix, numgroups)
      else
         call mc_write_prcfil(idmedium, int_mb(k_groupnums), prcfil,
     &     numgroups)
         call mc_build_prcfil(prcfil2, myid)
         write(fprefix2,'(256(a))') (' ', i=1,256)
         if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix2))
     &     call errquit('mc_main:rtdb get file_prefix failed'
     &      ,0,RTDB_ERR)
            write(mc_dir,'(256(a))') (' ', i=1,256)
         if (.not. rtdb_cget(rtdb, 'mc:directory',1,mc_dir))
     &      mc_dir(1:2)='./'
           write(fprefix,'(256(a))') (' ', i=1, 256)
           write(fprefix,'(3(a))')
     &        mc_dir(1:inp_strlen(mc_dir)),
     &        '/',
     &        fprefix2(1:inp_strlen(fprefix2))
      endif
      write(iw,*) 'DONE SETTING I/O files'
C-- the stdout is now redirected to other .out output files for each trajectory
C from here we can start mc_loop 
c      goto 100
      call task_energy(rtdb)
      if (ga_nodeid().eq.0) then
        iw=6   
        write(iw,*)  '\n\t\t METROPOLIS MONTE CARLO RECROSSING \n'
      endif
      
      write(iw,*) 'Subgroup is: ',util_sgroup_mygroup()
      dbug=.false. 
      forward=.true.
      backward=.false.
      write(iw,*)  "groupnum :"
      do i=1,numgroups
         write(iw,*)  int_mb(k_groupnums+i-1)
      enddo
       status=mc_data_set_ifirc(.true.)  
       status=mc_driver(rtdb,naccept)
      write(iw,*) 'Number of initial points generated: ', naccept  
c      go to 100 
      ircflag=.true.
      mdtemp10 =10.0 
      mdtemp = 298.15

      natom=mc_data_get_natom()

      call md_data_set_natom(natom)
      call md_data_allocate()
      call md_set(rtdb) ! sets md_temp, md_steps and mc_stepsize
      mdtemp=md_data_get_temp()
      tag=md_data_get_i_t()
      charge=md_data_get_i_q()
      c_start=md_data_get_i_c()
      vel=md_data_get_i_v()
      acc=md_data_get_i_a()
c  100 continue     
c --- not sure if I need a ga at all
      write(iw,*) '\nStarting the IRC trajectory\n'

c  ---- do the IRC trajectory first,  at low temp
C  ---- this requires quenching of velocities so that Ek is close to 0
      status=md_data_set_temp(mdtemp10)
C --- set vel and acc      

      status=mc_init(rtdb,grad,ircflag)

      c_array=mc_data_get_i_c_array()
      vel_init=mc_data_get_i_in_vel()
      acc_init=mc_data_get_i_in_acc()
      i_sc=mc_data_get_i_sc()
      i_s=mc_data_get_i_s()
      prp=mc_data_get_i_prp()
      mcsteps=mc_data_get_mcsteps()

      do k=1,3*natom
          dbl_mb(c_start+k-1)=dbl_mb(c_array+k-1)
          dbl_mb(vel+k-1)=dbl_mb(vel_init+k-1)
          dbl_mb(acc+k-1)=dbl_mb(acc_init+k-1)
      enddo

      if(dbug ) then
       write(iw,*) 'In MC main coordinates'
       do i=1,naccept+1
         write(iw,*) 'point ',i
         do k=1,3*natom
            write(iw,*) dbl_mb(c_array+(i-1)*3*natom+k-1)
         enddo
         write(iw,*) 'vel ',i
         do k=1,3*natom
            write(iw,*) dbl_mb(vel_init+(i-1)*3*natom+k-1)
         enddo
         write(iw,*) 'acc ',i
         do k=1,3*natom
            write(iw,*) dbl_mb(acc_init+(i-1)*3*natom+k-1)
         enddo

       enddo
      endif
C  --- run the IRC traj forward 

      status=mc_data_set_forward(.true.)
      status=mc_data_set_backward(.false.)
      status=mc_data_set_side_prev(0)
      status=mc_data_set_nxing(0)
      status=mc_data_set_trajnum(0)

      status=md_driver(rtdb,grad,thr)      
      status=mc_data_set_side_prev(0)
      status=mc_data_set_nxing(0)
C --- run the IRC traj backward
      call md_set(rtdb)
      do k=1,3*natom
          dbl_mb(c_start+k-1)=dbl_mb(c_array+k-1)
          dbl_mb(vel+k-1)=-dbl_mb(vel_init+k-1)
          dbl_mb(acc+k-1)=-dbl_mb(acc_init+k-1)
      enddo

      status=mc_data_set_forward(.false.)
      status=mc_data_set_backward(.true.)
      status=md_driver(rtdb,grad,thr) 

      call ga_sync()
C --- get geometry from the big array in mc_driver
C --- the first one is always the TS
      write(iw,*) 'Done with IRC'
      ircflag=.false.
      status=mc_data_set_ifirc(.false.)
C --- generate initial velocities for all trajectories
      status=mc_init(rtdb,grad,ircflag)

      do i=1,naccept+1
        write(iw,*) 'traj, prp: ',i,dbl_mb(prp+i-1)
      enddo

      write(iw,*) 'STARTING MD TRAJECTORIES '

      status=mc_data_set_forward(.true.)
      status=mc_data_set_backward(.false.)

C --- START LOOP OVER MD TRAJ - PARALLEL 
      call ga_sync()
      ncpu=numgroups

      do i=1,naccept+1
C ---- decide where to run the traj 
      if(util_sgroup_mygroup().eq.mod(i,ncpu)+1) then

C ---- set the appropriate coordinate             
        call md_set(rtdb)
        do k=1,3*natom
          dbl_mb(c_start+k-1)=dbl_mb(c_array+(i-1)*3*natom+k-1)
          dbl_mb(vel+k-1)=dbl_mb(vel_init+(i-1)*3*natom+k-1)
          dbl_mb(acc+k-1)=dbl_mb(acc_init+(i-1)*3*natom+k-1)
        enddo
        if (dbug ) then
          write(iw,*) 'Coordinates:'
           do k=1,3*natom
             write(iw,*) dbl_mb(c_start+k-1)
           enddo
          write(iw,*) 'Velocieties:'
           do k=1,3*natom
             write(iw,*) dbl_mb(vel+k-1)
           enddo
          write(iw,*) 'Acceleration:'
           do k=1,3*natom
             write(iw,*) dbl_mb(acc+k-1) 
           enddo
        endif
        do  j=1,3
          if (forward.and. .not.backward) then
            write(iw,*) 'Forward: traj',i  

            status=mc_data_set_side_prev(0)
            status=mc_data_set_nxing(0)
            status=mc_data_set_trajnum(i)
            status=md_driver(rtdb,grad,thr)
            forward=.false.
            backward=.true. 
            status=mc_data_set_forward(forward)
            status=mc_data_set_backward(backward)
            status=mc_data_set_side_prev(0)
            status=mc_data_set_nxing(0)

          elseif(backward.and. .not.forward) then
            write(iw,*) 'Backward: traj',i
C --- reset the geometry to c_start
C --- get the initial conditions generated in the forward run  
            call md_set(rtdb)
            do k=1,3*natom
              dbl_mb(c_start+k-1)=dbl_mb(c_array+(i-1)*3*natom+k-1)
              dbl_mb(vel+k-1)=-dbl_mb(vel_init+(i-1)*3*natom+k-1)
              dbl_mb(acc+k-1)=-dbl_mb(acc_init+(i-1)*3*natom+k-1)
            enddo
            status=mc_data_set_side_prev(0)
            status=mc_data_set_nxing(0)
            status=md_driver(rtdb,grad,thr)
            forward=.false.
            backward=.false.
            status=mc_data_set_side_prev(0)
            status=mc_data_set_nxing(0)

          write(iw,9000) 'traj ifirc forside0 backside0 forside 
     &                    backside', 
     &  i,
     &  mc_data_get_ifirc(),
     &  mc_data_get_forside0(),mc_data_get_backside0(),
     &  mc_data_get_forside(),mc_data_get_ backside()

            status=mc_data_set_forward(forward)
            status=mc_data_set_backward(backward)
          else
            forward=.true.
            backward=.false.
            status=mc_data_set_forward(forward)
            status=mc_data_set_backward(backward)

C --- frajectory finished, moved to the other trajectory  
C --- here you would also check if the trajectory was reactive or not 
        endif  
        enddo 
c      else
c         write(*,*) 'SKIP TRAJECTORY ',i,' node ',ga_nodeid()
      endif
c      write(*,*) 'LOOP ',i,' node ',ga_nodeid()
      enddo
c      write(*,*) ' Out of loop. Node id :',ga_nodeid()
C --- loop over trajectories finished
      call ga_sync()
      call ga_pgroup_set_default(idbig) 
      call ga_sync()
      call ga_igop(MT_F_INT,int_mb(i_sc),(mcsteps+2)*2,'+')
      call ga_sync()
      write(iw,*) ' Recrossing Results:'


      nonreact=0
      reactive=0
      reverse=0

      iw=6
      if (ga_nodeid().eq.0) then
       do i=1,naccept+1
        write(iw,9001) 'traj  sidef sideb ',i, int_mb(i_sc+2*i), 
     &               int_mb(i_sc+2*i+1) 
        if (int_mb(i_sc+2*i).eq.int_mb(i_sc+2*i+1)) then
           nonreact=nonreact+1
           dbl_mb(prp+i-1)=0.0 
           write(iw,*) 'trajectory number for non-reactive traj: ', i 
        else if (int_mb(i_sc+2*i).eq.mc_data_get_forside0()) then
           reactive=reactive+1
           dbl_mb(prp+i-1)=0.5*dbl_mb(prp+i-1)
        else
           reverse=reverse+1
           dbl_mb(prp+i-1)=-0.5*dbl_mb(prp+i-1)
        endif
       enddo

      j=0 
      avg=0 

       do i=0,mcsteps
        if (int_mb(i_s+i).eq.1) then
           j=j+1
        endif 
        write(iw,*) 'adding prp :',j,' val: ',dbl_mb(prp+j-1)    
        avg=avg+dbl_mb(prp+j-1)
       enddo

      avg=avg/(mcsteps+1)
      temp = mc_data_get_temp()
      rc =avg/sqrt((temp*boltz)/(twopi))

       write(iw,*) 'non-reactive, reactive,reverse traj: ',nonreact,
     &            reactive,reverse
       write(iw,*)  ' Average prp value: ',avg
       write(iw,*)  ' Recrossing coef: ',rc
       endif
c  100 continue 
      call ga_pgroup_set_default(idbig)
      call ga_sync()
      call mc_data_free_all()
      call md_data_free_all() 
c --- close files that were open
      close(iw)
      close(12)
      call util_sgend(rtdb)
      mc_main = .true.
 9000 format(A,I4,L3,4(I4))
 9001 format(A,3(I4))
      return    
      end

C***********************************************************************
      subroutine mc_setup_ga
     & (rtdb, myid, idbig, inodesbig,
     & idmedium, idzero, g_pr, groupnums, numgroups) 
C***********************************************************************
      Implicit none
C Include Files
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "msgids.fh"
#include "msgtypesf.h"
#include "util.fh"
#include "util_sgroup.fh"
#include "subgr.fh"
#include "pstat.fh"

C Variable Declarations
      integer rtdb ! input
C      integer nob, nspc
      integer myid, idbig, inodesbig, numgroups ! input
C Depreciated
      logical status ! internal use
C Depreciated
      integer idmedium
      integer ld(2)! internal use
      integer g_proc, procnums, groupnums, bigproc !p_proc is internal
C     procnums is internal
      integer i, j, idzero ! i and j are internal indexes
      integer ndim, dims(2), chunk(2), g_pr !n(dim)s and chunck are internal
      external util_sggo
C Additional variables
      integer groups_want, array_cpu(1), method
      integer dir
C End Additions

C Dimensions
      dimension procnums(numgroups)
      dimension groupnums(numgroups)
      dimension bigproc(inodesbig,1)

C timers common
C Create big GAs
      if (.not.ga_create(MT_DBL, inodesbig, 1, "proc list",
     &     1, -1, g_proc))
     &   call errquit('dntmc_setup_ga:g_proc create error', 0,
     &                 GA_ERR)
      call ga_fill(g_proc, -1)

C   Create Subgroups
      if (.not.rtdb_get(rtdb, 'subgroups_number', mt_int, 1,
     &   groups_want))groups_want=ga_nnodes()
C      write(*,*) 'groups_want =',groups_want
C groups_want is the number of subgroups, so initialy we want each CPU to be subgroup
C  that is all done in the input with set "subgroups_number"
     
C Setting Method = 1
c     1 -- use groups_want to generate equal sized groups (array_cpu ignored)
C     Simplest.
c     2 -- turn each SMP box into a group (array_cpu and groups_want ignored)
C     This uses GA to tell it about the cluster.
c     3 -- use array_cpu(groups_want) to define number of nodes per group
c     4 -- use array_cpu(groups_want+nnodes) to define which nodes per group
C     This is just option 3, but you get to lay the groups out exactly.
      method = 1
      array_cpu(1) = 0
      dir = 1 ! Write group rtdb's in scratch directories

      call util_sggo(rtdb,groups_want,method,array_cpu,dir)

C      if (util_print('debug',print_debug)) then
        write(6,*)'Now in Subgroups'
        call flush(6)
C      endif

C   I this is the handle to get into subgroups subgroup
      idmedium = ga_pgroup_get_default() ! Now in Subgroups
C   This is the handele to a group of 0th cpus
C    for example: 6 cpus in 3 groups, will have zerogroup(0,2,4)
      idzero = util_sgroup_zero_group()
C        write(6,*) 'idbig,idmedium, idzero:',idbig,idmedium,idzero

C Get Back to Big Group
      call ga_pgroup_set_default(idbig) ! Now in Big Group
C Create Processor Zero and Group Processor Lists
      if (ga_pgroup_nodeid(idmedium).eq.0) then
        ld(1)  = 1   !must be physical dimension of local array
        ld(2)  = 1
        call ga_put(g_proc,myid+1,myid+1,1,1,util_sgroup_mygroup(),ld)
        call ga_sync()
        ld(1)  = inodesbig  !must be physical dimension of local array
        ld(2)  = 1
        call ga_get(g_proc,1,inodesbig,1,1,bigproc,ld)
        j = 0
        do i = 1, inodesbig
C          write(*,*) 'big proc (i,1): ',bigproc(i,1)
          if (bigproc(i,1).ne.-1) then
            j = j + 1
            groupnums(j) = bigproc(i,1)
            procnums(j) = i-1
          endif
        enddo
        if (j.ne.util_sgroup_numgroups())
     &    call errquit('mc_setup_ga:zero node creation problem',
     &      j, UNKNOWN_ERR)
      else
        call ga_sync()
      endif
c       write(6,*)'Now broad:',idmedium,msg_dntmc3+MSGINT, procnums
c       write(6,*)'Now broad:',mitob(numgroups)
C Broadcase Results, array of groups( starting from 0)
      call ga_pgroup_brdcst(idmedium,msg_dntmc3+MSGINT, procnums,
     &     mitob(numgroups), 0)
C Broadcase Results, array of groups( starting from 1)
c       write(6,*)'Now broadcasting,again'
      call ga_pgroup_brdcst(idmedium,msg_dntmc7+MSGINT, groupnums,
     &     mitob(numgroups), 0)
c       write(6,*)'Now in destroying list'
C Destroy GA "proc list"
      if (.not. ga_destroy(g_proc)) call errquit('mc_setup_ga:
     &ga_destroy(g_proc) failed', GA_ERR)

C Create Zeros GAs
C this may not be needed, or will be needed for setting hared data for a subgroup

C      if (ga_pgroup_nodeid(idmedium).eq.0) then
C        call ga_pgroup_set_default(idzero)
C        write(*,*) 'In zero group nnodes ',ga_nnodes()
C          ndim = 2
C Changed to dimension (nob,nspc,numgroups)
C          dims(1) = NOB*nspc
C          dims(2) = util_sgroup_numgroups()
C          chunk(1) = dims(1)
C          chunk(2) = 0
C          g_pr = ga_create_handle()
C          call ga_set_data(g_pr,ndim,dims,MT_F_DBL)
C          call ga_set_chunk(g_pr,chunk)
C          status = ga_allocate(g_pr)
C          call ga_zero(g_pr)
C       endif
       write(6,*)'Now starting subgroups '
C Start Subgroups
       call ga_pgroup_set_default(idmedium) ! Now in Subgroups
c       write(*,*) 'In isubgroups nnodes ',ga_nnodes()
       call ga_sync()
C Test write
c       write(*,*) 'Depth: ',depth
c       write(6,*) '@ proc ',myid,' gr ',idmedium,'gr proc',ga_nodeid()
c       write(6,*) '@ proc ',myid,' gr2 ',my_ga_grp(depth),
c     & ' @ with zero group ',util_sgroup_zero_group()
c       if (myid .eq. 0) then
c           write(6,*)'@ group and proc arrays:numgroups',numgroups
c           do i=1,numgroups
c             write(6,*)'@ ',groupnums(i),procnums(i)
c           EndDo
c       endif

       return
       end


C***********************************************************************
      subroutine mc_setup_group
     &(rtdb, idmedium, myid, groupnums, prcfil, prcfil2, fprefix,
     & numgroups)
C gets the appropriate path for the output files and opens them
C***********************************************************************
      Implicit none

C Include Statements
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "subgr.fh"
#include "util.fh"
#include "util_sgroup.fh"

C Variable declarations
      integer rtdb ! input
      integer myid ! input
      integer idmedium ! input
      integer numgroups ! input
      integer groupnums(numgroups) ! input
      character*4 prcfil, prcfil2
      integer i
C Indect output of file units
      character*256 fprefix, fprefixcat ! internal use only !fprefix output
      character*256 fprefix2, mc_dir

C Begin Main Program
      call mc_write_prcfil(idmedium, groupnums, prcfil, numgroups)
      call mc_build_prcfil(prcfil2, myid)
      write(fprefix2,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix2))
     &     call errquit('dntmc_setup_group:rtdb get file_prefix failed'
     &     ,0,RTDB_ERR)
      write(mc_dir,'(256(a))') (' ', i=1,256)
      if (.not. rtdb_cget(rtdb, 'mc:directory',1,mc_dir))
     &   mc_dir(1:2)='./'
      write(fprefix,'(256(a))') (' ', i=1, 256)
      write(fprefix,'(3(a))')
     &     mc_dir(1:inp_strlen(mc_dir)),
     &     '/',
     &     fprefix2(1:inp_strlen(fprefix2))

C Open file Units
C general output ! Only group Zeros write
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(3(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.xyz.',
     &      prcfil(1:4)

c -- open a file for trajectory for each subgroup
c      if(util_print('MCout',print_low)) then
      OPEN(UNIT = 12,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')
c      endif

      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(3(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.out.',
     &      prcfil(1:4)

C  ---- added this so that the output file is closed and output is open for each 
C  ---- subgroup
c      if(util_print('out',print_low)) then
      OPEN(UNIT = 6,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')
c      endif


      return
      end

C***********************************************************************
      subroutine mc_write_prcfil(idmedium,groupnums,prcfil,numgroups)
C finds the appropriate groupid, to get the right output file
C***********************************************************************
      Implicit none

C Include Statements
#include "subgr.fh"
#include "errquit.fh"
#include "util_sgroup.fh"

C Variable Declarations
      integer i
      integer numgroups
      integer idmedium, groupnums(numgroups)
      integer groupid
      character*4 prcfil ! only output

C Main Program
      groupid = -1
      do i=1, util_sgroup_numgroups()
        if (groupnums(i) .eq. util_sgroup_mygroup()) then
          groupid = i
        endif
      enddo

      if (groupid .eq. -1)
     & call errquit('dntmc_write_prcfil:failed to allocate idgroup',0,
     &   GA_ERR)

       call mc_build_prcfil(prcfil, groupid)

C End Main Program
      return
      end

***********************************************************************
      subroutine mc_build_prcfil(prcfil, i)
C formats the name of the procesor file: ex. 0001,0010,0100
C so you can have up to 9999 cpus and files
C***********************************************************************
      Implicit none

C  Variable Declarations
      integer i
      character*4 prcfil ! only output

C  Main Program
      write(prcfil(1:4), '(i4)') i
      if (i .le. 9) then
         prcfil(1:3) = '000'
      endif
      if (i .le. 99) then
         prcfil(1:2) = '00'
      endif
      if (i .le. 999) then
         prcfil(1:1) = '0'
      endif

C  End Main Program
      return
      end

