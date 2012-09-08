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
#include "util_sgroup.fh"  
#include "const_data.fh"
#include "eaf.fh"
#include "subgr.fh"
      logical mc_main
      logical md_driver,mc_driver,mc_init
      integer mcsteps
      logical grad, task_energy
      external grad, task_energy

      logical status
      integer natom 
      integer c_start
      integer vel_init,acc_init
      logical forward,backward   ! tells you which part of the traj is being called
      logical mc_data_set_ifirc, mc_data_get_ifirc
      integer mc_data_get_forside, mc_data_get_forside0
      integer mc_data_get_backside, mc_data_get_backside0
      integer i,j,k,naccept 
      integer c_array, vel, acc, prp, prp_tmp,sidef,sideb 
      integer mc_data_get_i_c_array, md_data_get_i_v, md_data_get_i_a 
      integer md_data_get_i_c, mc_data_get_i_in_vel,mc_data_get_i_in_acc
      integer mc_data_get_natom, mc_data_get_i_prp
      integer  mc_data_get_i_prp_tmp
      integer tag,charge,md_data_get_i_t,md_data_get_i_q
      character*32 thr
      logical ircflag,mc_data_set_forward,mc_data_set_backward
      integer i_sc, mc_data_get_i_sc, nonreact,reactive,reverse
      integer i_s, mc_data_get_i_s,mc_data_get_mcsteps
      double precision avg,rc,rc1,temp,mc_data_get_temp
      logical md_data_set_temp
      logical mc_data_set_side_prev, mc_data_set_nxing
      logical mc_data_set_trajnum
      double precision mdtemp,mdtemp10, md_data_get_temp
      integer md_data_get_nsteps, nsteps
      logical dbug
c<<<<<<< .mine
      logical restart_dump

c=======
      character*(nw_max_path_len) fnameout
      integer tunita
c
c>>>>>>> .r18494
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
      mc_main=.false.
      iw=6
      call ga_sync()
C     INFO ON PROCESSES in the big group
      myid = ga_nodeid()
      idbig = ga_pgroup_get_default()
      inodesbig = ga_nnodes() 

      if (myid.eq.0) then
          open(unit=ir,file='restart.pun',status='replace')
      endif
      if (.not. rtdb_get(rtdb,'subgroups_number', MT_INT, 1, numgroups))
     &    numgroups = ga_nnodes() 

C ***** Can't have more subgroups then procs *****
      if (numgroups.gt.ga_nnodes()) then
          numgroups = ga_nnodes()
          write (iw,*) 'Requested number of subgroups too big, 
     &    reseting to:', numgroups
      endif
C     Allocate memory for array of group numbersi, only if more then 1 subgroup requested
      if (numgroups.gt.1) then
      if (.not. ma_push_get(MT_INT, numgroups, 'mc: groupnums',
     &    l_groupnums, k_groupnums)) 
     &    call errquit('mc_allocate_arrays: error groupnums',0,MA_ERR)

     
      if (ga_nodeid().eq.0) then
          write(iw,*) 'SPLITTING INTO SUBGROUPS'
      endif
      call ga_sync()
      call util_flush(iw)
      close(iw) 
C    SPLITTING INTO SUBGROUPS

C   Close the output file and reopen it as series of ".out" files
      call mc_setup_ga
     &     (rtdb, myid, idbig, inodesbig,
     &     idmedium, idzero, g_pr, int_mb(k_groupnums),numgroups) ! Now in Subgroups

C Initialize group related information
C mostly file opening for each group,etc. 
      endif
      call mc_setup_group(numgroups)  !deals with file i/o
c     &     (rtdb, myid, int_mb(k_groupnums),prcfil, prcfil2,
c     &     fprefix, numgroups)

      write(iw,*) 'DONE SETTING I/O files for subgroups'
C-- the stdout is now redirected to other .out output files for each trajectory
C from here we can start mc_loop 
c      goto 100
c      status = task_energy(rtdb)
c      goto 100 

      write(iw,*) 'testing energy calculation\n'
      status = task_energy(rtdb)

      write(iw,*)  '\n\t\t METROPOLIS MONTE CARLO RECROSSING \n'
      
      write(iw,*) 'Subgroup is: ',util_sgroup_mygroup()
      dbug=.false. 
      forward=.true.
      backward=.false.

C********* MC driver ******

C******* RESTART point 1******
      if(.not.mc_driver(rtdb,naccept)) call errquit('mc_driver:failed',
     &    UNKNOWN_ERR)
      write(iw,*) 'Number of initial points generated: ', naccept  

      if(util_sgroup_mygroup().eq.1.and.ga_nodeid().eq.0) then
         status=restart_dump()
      endif
C******** End of MC driver *******
c      go to 100 
      ircflag=.true.
      status=mc_data_set_ifirc(.true.)
      mdtemp10 =10.0 
      mdtemp = 298.15

      natom=mc_data_get_natom()

      call md_data_set_natom(natom)
      call md_data_allocate()
      call md_set(rtdb) ! sets md_temp, nsteps and mc_stepsize
      mdtemp=md_data_get_temp()
      tag=md_data_get_i_t()
      charge=md_data_get_i_q()
      c_start=md_data_get_i_c()
      vel=md_data_get_i_v()
      acc=md_data_get_i_a()
c  100 continue     
      write(iw,*) '\nStarting the IRC trajectory\n'

c  ---- do the IRC trajectory first,  at low temp
C  ---- this requires quenching of velocities so that Ek is close to 0
C  ---- quenching not implemented yet
      status=md_data_set_temp(mdtemp10)
C --- set vel and acc      

      status=mc_init(rtdb,grad,ircflag)

      c_array=mc_data_get_i_c_array()
      vel_init=mc_data_get_i_in_vel()
      acc_init=mc_data_get_i_in_acc()
      i_sc=mc_data_get_i_sc()
      i_s=mc_data_get_i_s()
      prp=mc_data_get_i_prp()
      prp_tmp = mc_data_get_i_prp_tmp()
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

      if(.not.md_driver(rtdb,grad,thr)) 
     &  call errquit('Error in MDDriver',0,0)      
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

      if(.not.md_driver(rtdb,grad,thr))
     &  call errquit('Error in MDDriver',0,0)  

      call ga_sync()

C ****** RESTART point 2 ******* 
c ---- dump the data for restart
      if(util_sgroup_mygroup().eq.1.and.ga_nodeid().eq.0) then
         status=restart_dump()
      endif

      write(iw,*) 'Done with IRC'
      ircflag=.false.
      status=mc_data_set_ifirc(.false.)
C --- generate initial velocities for all trajectories
      status=mc_init(rtdb,grad,ircflag)

      do i=1,naccept+1
        write(iw,*) 'traj, prp: ',i,dbl_mb(prp+i-1)
        dbl_mb(prp_tmp+i-1) = dbl_mb(prp+i-1)
      enddo

      write(iw,*) 'STARTING MD TRAJECTORIES '

      status=mc_data_set_forward(.true.)
      status=mc_data_set_backward(.false.)

C --- START LOOP OVER MD TRAJ - PARALLEL 
      call ga_sync()
      ncpu=numgroups

C ****** RESTART point 3 ******* 
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
            
            if(.not.md_driver(rtdb,grad,thr))
     &       call errquit('Error in MDDriver',0,0)  
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

            if(.not.md_driver(rtdb,grad,thr))
     &         call errquit('Error in MDDriver',0,0)  
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
c      call ga_sync()
c      idmedium = ga_pgroup_get_default()     
c      write(iw,*) 'big,med,node',idbig,idmedium,ga_nodeid() 
c      flush(iw)
c      call ga_pgroup_set_default(idbig)
* CCCCCC igop needs to sync all procs, so it makes the subgroups wait for eachother
* CCCCCc may need to have it as a global array        
c      call ga_igop(MT_F_INT,int_mb(i_sc),(mcsteps+1)*2,'max')
c      call ga_pgroup_set_default(idmedium)
c      if(util_sgroup_mygroup().eq.1.and.ga_nodeid().eq.0) then
c           status=restart_dump()
c      endif
      enddo
c      write(*,*) ' Out of loop. Node id :',ga_nodeid()
C --- loop over trajectories finished
      call ga_sync()
      call ga_pgroup_set_default(idbig) 
      call ga_sync()
      nsteps=md_data_get_nsteps()
      call ga_igop(MT_F_INT,int_mb(i_sc),(mcsteps+1)*2*nsteps,'max')
      call ga_sync()
      status=restart_dump()
      write(iw,*) ' Recrossing Results:'



      if (ga_nodeid().eq.0) then
c loop over each mdstep to get kappa for each step      
       rc1=0.0
       do k=0,nsteps-1
          nonreact=0
          reactive=0
          reverse=0
          do i=1,naccept+1
c      int_mb(i_sc+(trajnum-1)*2*nsteps+mdstep-1) = side
             sidef = int_mb(i_sc+2*(i-1)*nsteps+k)
             sideb = int_mb(i_sc+2*(i-1)*nsteps+nsteps+k)
             if(dbug)
     &         write(iw,9001) 'step traj sidef sideb ',k,i,sidef,sideb
             if (sidef.eq.sideb) then
                nonreact=nonreact+1
                dbl_mb(prp_tmp+i-1)=0.0 
c                if (dbug)
                 write(iw,*) 'tr number for non-reactive traj: ', i 
             else if (sidef.eq.mc_data_get_forside0()) then
                dbl_mb(prp_tmp+i-1)=0.5*dbl_mb(prp_tmp+i-1)     
c                 reactive=reactive+1
             else
c                reverse=reverse+1
                dbl_mb(prp_tmp+i-1)=-0.5*dbl_mb(prp_tmp+i-1)
             endif
             if (dbl_mb(prp_tmp+i-1).gt.0.0) then
                reactive=reactive+1
             else if (dbl_mb(prp_tmp+i-1).lt.0.0) then
                reverse=reverse+1
             endif
          enddo

          j=0 
          avg=0 


          do i=0,mcsteps
            if (int_mb(i_s+i).eq.1) then
              j=j+1
            endif
c            if(dbug) then 
              write(iw,*) 'Mc step: ',i,'adding: ',j,' val: ',
     &        dbl_mb(prp_tmp+j-1)
c            endif
              avg=avg+dbl_mb(prp_tmp+j-1)
          enddo

          avg=avg/(mcsteps+1)
          temp = mc_data_get_temp()
          rc=avg
c          rc =avg/sqrt((temp*boltz)/(twopi))

          if(k.eq.0) rc1=rc

          write(iw,9001) 'non-reactive, reactive,reverse traj: ',
     &         nonreact,reactive,reverse
          write(iw,*)  ' Average prp value: ',avg
          write(iw,*)  'Analytical avg: ',sqrt((temp*boltz)/(twopi)) 
          write(iw,9002)  'MD Step,  Recrossing coef: ',k+1,rc/rc1
          write(iw,9002)  'MD Step,  Recrossing coef anl: ',k+1,
     $          rc/sqrt((temp*boltz)/(twopi))          
C   --- reset the original prp values         
          do i=1,naccept+1
             dbl_mb(prp_tmp+i-1) = dbl_mb(prp+i-1)
          enddo
        enddo  ! end of loop over nsteps
       endif
c  100 continue 
      call ga_pgroup_set_default(idbig)
      call ga_sync()
      call mc_data_free_all()
      call md_data_free_all() 
c --- close files that were open
      close(iw)
      close(12)
      close(ir)
      call util_sgend(rtdb)
      mc_main = .true.
 9000 format(A,I6,L3,4(I6))
 9001 format(A,3(I6))
 9002 format(A,I6,F10.6)
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

      if (util_print('debug',print_debug)) then
        write(0,*)'Now in Subgroups'
        call flush(6)
      endif

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
     &ga_destroy(g_proc) failed', 0, GA_ERR)

C Create Zeros GAs
C this may not be needed, or will be needed for setting hared data for a subgroup

c       write(6,*)'Now starting subgroups '
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
      subroutine mc_setup_group(numgroups)
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
#include "global.fh"
C Variable declarations
c      integer rtdb ! input
c      integer myid ! input
c      integer numgroups ! input
c      integer groupnums(numgroups) ! input
c      character*4 prcfil, prcfil2
c      integer i
C Indect output of file units
c      character*256 fprefix, fprefixcat ! internal use only !fprefix output
c      character*256 fprefix2, mc_dir
      integer numgroups
      character*256 xyzname, outname
      character*256 pxyzname, poutname
     
C Begin Main Program
c      call mc_write_prcfil(groupnums, prcfil, numgroups)

c     call mc_build_prcfil(prcfil2, myid)
c      write(fprefix2,'(256(a))') (' ', i=1,256)
c      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix2))
c     &     call errquit('dntmc_setup_group:rtdb get file_prefix failed'
c     &     ,0,RTDB_ERR)
c      write(mc_dir,'(256(a))') (' ', i=1,256)
c      if (.not. rtdb_cget(rtdb, 'mc:directory',1,mc_dir))
c     &   mc_dir(1:2)='./'
c      write(fprefix,'(256(a))') (' ', i=1, 256)
c      write(fprefix,'(3(a))')
c     &     mc_dir(1:inp_strlen(mc_dir)),
c     &     '/',
c     &     fprefix2(1:inp_strlen(fprefix2))
C Open file Units
C general output ! Only group Zeros write
      if (ga_nodeid().eq.0 ) then 
          call util_file_name('xyz', .false., .false., xyzname)
          call util_file_name('out', .false., .false., outname)
          call util_pname(xyzname,pxyzname)
          call util_pname(outname,poutname) 
          write(6,*) 'opening files', pxyzname, poutname 
c          write(fprefixcat,'(256(a))') (' ', i=1, 256)
c          write(fprefixcat,'(3(a))')
c     &      fprefix(1:inp_strlen(fprefix)),
c     &      '.xyz.',prcfil(1:4)
c           write(6,*) 'Filename', xyzname
c -- open a file for trajectory for each subgroup
           if(numgroups.gt.1) then 
             OPEN(UNIT = 12,FILE = pxyzname,STATUS = 'REPLACE')

C  ---- added this so that the output file is closed and output is open for each 
C  ---- subgroup
             OPEN(UNIT = 6,FILE = poutname,STATUS = 'REPLACE')
           else
               OPEN(UNIT = 12,FILE = xyzname,STATUS = 'REPLACE')
               OPEN(UNIT = 6,FILE = outname,STATUS = 'REPLACE')
           endif
      else
C all non-zeros write to dev/null
            OPEN(6,FILE='/dev/null',STATUS='UNKNOWN')
            OPEN(12,FILE='/dev/null',STATUS='UNKNOWN')
      endif 


      return
      end

C***********************************************************************
      subroutine mc_write_prcfil(groupnums,prcfil,numgroups)
C finds the appropriate groupid, to get the right output file
C not using idmedium, need to remove it
C***********************************************************************
      Implicit none

C Include Statements
#include "subgr.fh"
#include "errquit.fh"
#include "util_sgroup.fh"

C Variable Declarations
      integer i
      integer numgroups
      integer groupnums(numgroups)
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

      function restart_dump() 
      Implicit none
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
#include "util_sgroup.fh"  
#include "const_data.fh"
#include "subgr.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh" 
      logical restart_dump
C: save the array of initial velocities and initial structures and initial prp values
C: save the array with 1 and 0 for iniital values
C: before the run, save values to -1, so at restart we know which ones to get.      
C: restart file has handle ir=11 (set as parameter in const_data.fh
      integer prp, vel, acc, coor, sidecnt, mcsteps, naccept, natom
      integer nsteps
      integer mc_data_get_i_prp !property (naccept+2)
      integer mc_data_get_i_sc  !side count (mcsteps+1)*2*nsteps
      integer mc_data_get_i_in_vel ! initial velocities (naccept+2)*3*natom
      integer mc_data_get_i_in_acc ! initial accelerations  (naccept+2)*3*natom
      integer mc_data_get_i_c_array ! initial coordinates (mcsteps+1)*3*natom
      integer mc_data_get_natom     ! number of atoms
      integer mc_data_get_mcsteps   ! mc steps mcsteps+1
      integer mc_data_get_naccept   ! number of accepted steps
      integer myid, mysubgr
      integer md_data_get_nsteps    ! number of steps in MD trajectory
      integer i
C*******  double check the size of the array for vel,acc and prp****
      restart_dump=.false.
      natom=mc_data_get_natom()
      mcsteps=mc_data_get_mcsteps()
      naccept=mc_data_get_naccept()
      prp=mc_data_get_i_prp()
      vel=mc_data_get_i_in_vel()
      acc=mc_data_get_i_in_acc() 
      coor=mc_data_get_i_c_array()
      sidecnt=mc_data_get_i_sc()  
      nsteps = md_data_get_nsteps()
C:   Check again if it is subgroup 1, node 0 before you write
      myid = ga_nodeid()
      mysubgr = util_sgroup_mygroup() 
c      write(ir,*) 'writing to restart file'
      if (myid.eq.0.and.mysubgr.eq.1) then
          rewind(ir)
          write(ir,*) 'natom,mcsteps,naccept,nsteps'
          write(ir,*) natom,mcsteps,naccept,nsteps    
C figure out a better print      
          write(ir,*) 'structures'
          call print_array_dbl(dbl_mb(coor),(mcsteps+1)*3*natom,ir)
          write(ir,*) 'sidecount' 
          call print_array_int(int_mb(sidecnt), (mcsteps+1)*2*nsteps,ir)
          write(ir,*) 'prp'
          call print_array_dbl(dbl_mb(prp),(naccept+1),ir)
          write(ir,*) 'vel'
          call print_array_dbl(dbl_mb(vel),(naccept+1)*3*natom,ir)
          write(ir,*) 'acc'
          call print_array_dbl(dbl_mb(acc),(naccept+1)*3*natom,ir)
          call util_flush(ir)
      endif 
      restart_dump=.true.
      return 
      end
C ******** End of writing restart file *********
C ******** Start reading restart file *********
      function restart_read()
      Implicit none
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
#include "util_sgroup.fh"  
#include "const_data.fh"
#include "subgr.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      logical restart_read
C: save the array of initial velocities and initial structures and initial prp values
C: save the array with 1 and 0 for iniital values
C: before the run, save values to -1, so at restart we know which ones to get.      
C: restart file has handle ir=11 (set as parameter in const_data.fh
      integer prp, vel, acc, coor, sidecnt, mcsteps, naccept, natom
      integer nsteps
      integer mc_data_get_i_prp !property (naccept+2)
      integer mc_data_get_i_sc  !side count (mcsteps+1)*2*nsteps
      integer mc_data_get_i_in_vel ! initial velocities (naccept+2)*3*natom
      integer mc_data_get_i_in_acc ! initial accelerations  (naccept+2)*3*natom
      integer mc_data_get_i_c_array ! initial coordinates (mcsteps+1)*3*natom
      logical mc_data_set_naccept   ! number of accepted steps
      logical md_data_set_nsteps
      integer myid, mysubgr
      integer i
      logical status
C*******  double check the size of the array for vel,acc and prp****
      restart_read=.false.
      prp=mc_data_get_i_prp()
      vel=mc_data_get_i_in_vel()
      acc=mc_data_get_i_in_acc()
      coor=mc_data_get_i_c_array()
      sidecnt=mc_data_get_i_sc()
C:   Check again if it is subgroup 1, node 0 before you write
c      write(ir,*) 'writing to restart file'
      rewind(ir)
      read(ir,*) natom,mcsteps,naccept,nsteps
C figure out a better print      
          call read_array_dbl(dbl_mb(coor),(mcsteps+1)*3*natom,ir)
          call read_array_int(int_mb(sidecnt), (mcsteps+1)*2*nsteps,ir)
          call read_array_dbl(dbl_mb(prp),(naccept+1),ir)
          call read_array_dbl(dbl_mb(vel),(naccept+1)*3*natom,ir)
          call read_array_dbl(dbl_mb(acc),(naccept+1)*3*natom,ir)
      status=mc_data_set_naccept(naccept)
      restart_read=.true.
      return
      end
C ******** End reading restart file **********



      subroutine print_array_int(array,asize,fileout)
      implicit none
      integer asize,fileout
      integer array(asize)
      write(fileout,*) array
      return
      end
      subroutine print_array_dbl(array,asize,fileout)
      implicit none
      integer asize,fileout
      double precision array(asize)
      write(fileout,*) array
      return
      end
      subroutine read_array_int(array,asize,filein)
      implicit none
      integer asize,filein
      integer array(asize)
      read(filein,*) array
      return
      end
      subroutine read_array_dbl(array,asize,filein)
      implicit none
      integer asize,filein
      double precision array(asize)
      read(filein,*) array
      return
      end


c $Id$
