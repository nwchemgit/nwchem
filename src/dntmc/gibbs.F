C$Id: gibbs.F,v 1.16 2007-08-03 23:01:19 d3p852 Exp $
C23456789012345678901234567890123456789012345678901234567890123456789012
C     modified on oct 3, 1999 for restart on mpp1
C     modified on oct 25, 1999 for random number generator restart
C     modified on may 10, 2004 for colony2a  
C     modified on april 15, 2005 for ran3

      logical function dntmc_GIBBS(rtdb)

      IMPLICIT none
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "msgtypesf.h"
C  added includes
#include "inp.fh"
C  End added includes

c
      integer rtdb  ! integer

      logical task_energy
      logical util_sgstart
      external util_sgstart
      external task_energy
      external util_get_rtdb_state

      double precision M,MTOT
      double precision r, rcm, rs
      double precision ri
      double precision rmin, rsim, rmax
      double precision e, k, esqrd, kb, pi
      double precision qse
      double precision qsqrd, sgmsq, epsln
      integer lndx, lndxx
      double precision lambda
      INTEGER NATOM
      double precision disp, alpha
      INTEGER ROTF
      INTEGER SPC
      INTEGER USPC
      integer l
      double precision abndx
      INTEGER AB
      double precision usum, usqrdsum
      double precision rmsdsum, rmsdavr 
      double precision rmsdsumsp, rmsdavrsp
      integer number
      double precision NUMBERSUM, NUMBERAVR
      double precision PHI, PHISUM, PHIAVR
      double precision rclcm
      INTEGER NMIN,NCURR,NMAX
      integer idum
      double precision seed
      double precision rnnum
      double precision rb, testme
      integer bigproc, procnums
      integer ngroups
#ifdef USE_SUBGROUPS
      integer util_sgroup_numgroups
      external util_sgroup_numgroups
#endif
      double precision NDA,NDSIG

      double precision temp, uu, uavr, usqrd, rsimsum, rsimavr
      double precision rsh, rsimsq, rmaxsq, rx, ry, rz, dist
      double precision rinit
      double precision sgm, aratiott, aratiotr, aratiotv
      double precision aratiotpr, aratiotpi, cv, rinc
      double precision percent, rest
      INTEGER NSPC, NDIM,g_pr,g_proc,DIMS,CHUNK,IICOUNT 
      INTEGER CONT1,CONT2
      integer myid
      integer mxnode, ld, natot, nob, iam
      integer iacceptt, iacceptr, iacceptv, iacceptpr, iacceptpi
      integer nmols, mcstps, lmin, lmax, mprnt, i, j, kp, node
      integer ir, ii, jj, iname, ilo, ihi, jlo, jhi
      integer idbig, inodesbig, idmedium, idzero

c  Variables specific to added restart
      integer procrestart
      integer rsmcstps
      integer lnow
      logical os, oc, or
      character*256 fprefix
      character*256 fprefixcat
      integer groupnums
      dimension groupnums(2000)
      character*4 prcfil2
c  Variables for save of old nda
      double precision nda_old
      dimension nda_old(100)
C  Variables for global sync
      integer g_stop, stoplist
      dimension stoplist(2000)
      logical trigger
C Added array for configuration files
      double precision DISTcalc(MXMOLS)
      integer molmax
C Variables for rweight distribution
      logical rweight 
C  End Variables

      character*4 prcfil
      CHARACTER*4 occtag
      CHARACTER*8 PREFIX
      CHARACTER*11 PREFIXX
      CHARACTER*11 PREFIX1
      CHARACTER*4 SUFFIX
      CHARACTER*4 SUFFIXX
      CHARACTER*4 SUFFIX1
      LOGICAL ACCEPTT
      LOGICAL ACCEPTR
      LOGICAL ACCEPTV
      LOGICAL ACCEPTPR
      LOGICAL ACCEPTPI
      LOGICAL STATUS

      PARAMETER (MXNODE = 2000)
      PARAMETER (PREFIX = "h2so4h2o")
      PARAMETER (SUFFIX = ".mol")
      PARAMETER (PREFIXX = "h2so4h2opdb")
      PARAMETER (SUFFIXX = ".inp")
      PARAMETER (PREFIX1 = "h2so4h2ocnt")
      PARAMETER (SUFFIX1 = ".out")

      COMMON/MASSES/M(MXSP,MXATOM),MTOT(MXSP)
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/COORDI/RI(MXMOLS,MXATOM,3)
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/INDEX/QSE(MXSP,MXATOM,3)
      COMMON/INDEXX/QSQRD(MXSP,MXSP,MXATOM,MXATOM),
     +SGMSQ(MXSP,MXSP,MXATOM,MXATOM),
     +EPSLN(MXSP,MXSP,MXATOM,MXATOM)
      COMMON/LAMBNDX/LNDX(MXMOLS),LAMBDA,LNDXX(9,9)
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/TRPS/DISP,ALPHA,ROTF
      COMMON/SPECIES/SPC(MXMOLS)
      COMMON/USPECIES/USPC(MXSP)
      COMMON/MSTP/L
      COMMON/ABENSBL/ABNDX(2)
      COMMON/AORB/AB
      COMMON/AVR/USUM,USQRDSUM
      COMMON/RMSDCOORD/RMSDSUM(MXMOLS),RMSDAVR(MXMOLS)
      COMMON/RMSDCOORDSP/RMSDSUMSP(MXSP),RMSDAVRSP(MXSP)
      COMMON/NUMBER1/NUMBER(100)
      COMMON/NUMBER2/NUMBERSUM(100),NUMBERAVR(100)
      COMMON/ANGLES/PHI(100),PHISUM(100),PHIAVR(100)
      COMMON/COORD2/RCLCM(3)
      COMMON/GRAND/NMIN,NCURR,NMAX
      COMMON/SEED/IDUM  
      COMMON/SEEDS/SEED(1)
      COMMON/RNNUMS/RNNUM(1)   
 
C      DIMENSION TESTME(100,MXNODE)
      COMMON/DIST1/RB(100),TESTME(100,2000)
      COMMON/DIST2/NDA(100),NDSIG(100)

      DIMENSION DIMS(2),CHUNK(2)
      DIMENSION ld(2)
      DIMENSION RINIT(MXMOLS,MXATOM,3)
      dimension bigproc(2000,1)
      dimension procnums(2000)

      SAVE NATOT

      dntmc_gibbs = .false.
      NOB = 100
      if (.not. rtdb_get(rtdb, 'dntmc:procrestart', MT_INT, 1,
     &    procrestart))
     &   procrestart = 0

      if (.not. rtdb_get(rtdb, 'dntmc:rweight', MT_LOG, 1, rweight))
     &   rweight = .false.

      if (procrestart .ne. 0) then
            write(fprefix,'(256(a))') (' ', i=1,256)
            if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix))
     &      call errquit('rtdb get file_prefix failed',0,RTDB_ERR)
            call dntmc_getspecinfo(rtdb, nmols, natom, rinit, uspc, nspc
     &      ,m)
            do i=1, nmols
              spc(i) = uspc(1) !assumes only one species
            enddo
        if (ga_nodeid() .eq. 0) then
          do i=1, procrestart
            call dntmc_build_prcfil(prcfil, i)
          call dntmc_copy_filertdb(prcfil, rtdb, nmols, fprefix)
          enddo
            if(.not. rtdb_delete(rtdb, 'dntmc:procrestart'))
     &        call errquit('rtdb failed to delete dntmc:procrestart'
     &        ,0,RTDB_ERR)
        endif
        call dntmc_sync(.false.)
        goto 9945
      endif

c
c  get the global node id number and group handle
c  this is before any of the subgroups are created
c
      myid = ga_nodeid()
      idbig = ga_pgroup_get_default()
c     idbig = ga_pgroup_get_world()
      inodesbig = ga_nnodes()

c     WRITE(6,*)'i made it into the code =',myid
c     CALL FLUSH(6)

C********************************************************
c
c  this is where the creation of subgroups probably 
c  needs to come in, each of the 0 nodes needs to be 
c  involved in the creation of the g_pr array
c
c     if (inodesbig.gt.1) then
      if (.not.ga_create(MT_DBL, inodesbig, 1, "proc list",
     &     1, -1, g_proc))
     &   call errquit('dntmc_gibbs:g_proc create error', 0,
     &                 GA_ERR)
      call ga_fill(g_proc, -1)
C  Global sync mutex creation
      if (.not. ga_create_mutexes(1))
     &   call errquit('dntmc_gibbs:mutex create error', 0, GA_ERR)
C  End Global sync mutex creation
       if (.not.util_sgstart(rtdb))
     &   call errquit('dntmc_gibbs:error starting subgroups',
     &                  0, UNKNOWN_ERR)
       idmedium = ga_pgroup_get_default()
#ifdef USE_SUBGROUPS
       ngroups = util_sgroup_numgroups()
#else
       ngroups = 1
#endif
c
c This next part is to set up the group that is only the node 0s
c in the medium groups
c

       call ga_pgroup_set_default(idbig)
       if (ga_pgroup_nodeid(idmedium).eq.0) then
c
c Only node 0s get in here and then they go to the world
c group to form the node 0 group.
c
          write(6,*) 'inside the ga_nodeid eq 0 clause'
          call flush(6)
C          call ga_pgroup_set_default(idbig)
C          write(6,*) 'after the set_default idbig'
C          call flush(6)
c
c create the proclist of only the CPU 0s
c
          ld(1)  = 1   !must be physical dimension of local array
          ld(2)  = 1
          call ga_put(g_proc,myid+1,myid+1,1,1,idmedium,ld)
C         call ga_sync()
          call dntmc_sync(.false.)

c
c all 0's get the full list
c
          ld(1)  = 2000   !must be physical dimension of local array
          ld(2)  = 1
          call ga_get(g_proc,1,inodesbig,1,1,bigproc,ld)
          call dntmc_sync(.false.)
          write(6,*) 'after ga_get of the node ids'
          call flush(6)
          j = 0
          do i = 1, inodesbig
             if (bigproc(i,1).ne.-1) then
                j = j + 1
                groupnums(j) = bigproc(i,1)
                procnums(j) = i-1
             endif
          enddo
c
c check to make sure that we got the correct number of processes
c
          write(6,*) 'j and ngroups ',j,ngroups,myid
          if (j.ne.ngroups)
     &      call errquit('dntmc_gibbs:zero node creation problem',
     &        j, UNKNOWN_ERR)
       else
c
c All nodes have to be back at the large group before the
c 0 node group can be formed.
c
          write(6,*) 'inside the ga_nodeid eq 0 else clause'
          call flush(6)
C          call ga_pgroup_set_default(idbig)
C          write(6,*) 'after the set_deafault idbig'
C          call flush(6)
C          call ga_sync()
          call dntmc_sync(.false.)
          call dntmc_sync(.false.)
       endif 
c
c Everyone need to get the processors involved in idzero and do the create
c

C       call ga_sync()
C       call dntmc_sync(.false.)
C       call ga_brdcst(msg_dntmc3+MSGINT, procnums, mitob(ngroups), 0)
C       call ga_brdcst(msg_dntmc7+MSGINT, groupnums, mitob(ngroups), 0)
       call ga_pgroup_brdcst(idmedium,msg_dntmc3+MSGINT, procnums,
     &                       mitob(ngroups), 0)
       call ga_pgroup_brdcst(idmedium,msg_dntmc7+MSGINT, groupnums, 
     &                       mitob(ngroups), 0)
C       call ga_sync()
       call dntmc_sync(.false.)

       write(6,*) 'before the ga_pgroup_create idzero call',
     &    procnums(1), procnums(2), ngroups, ga_nodeid(), myid
       idzero = ga_pgroup_create(procnums, ngroups)
       write(6,*) 'after the ga_pgroup_create idzero call'
c
c Now set the default group to idmedium so that we can get to the node
c zeros again!
c
C       call ga_pgroup_set_default(idmedium)
C        call ga_pgroup_set_default(idzero)
       if (ga_pgroup_nodeid(idmedium).eq.0) then
           call ga_pgroup_set_default(idzero)
C          write(6,*) 'after the ga_pgroup_set_default idzero call'

          NOB = 100  !number of bins

          ndim = 2

          dims(1) = NOB
          dims(2) = ngroups
c           dims(2) = ga_nnodes()

          chunk(1) = dims(1) !if = dims(1), gives distribution by vertical strips
          chunk(2) = 0       !if < 1, causes dimension to be distributed evenly


          g_pr = ga_create_handle()
          call ga_set_data(g_pr,ndim,dims,MT_F_DBL)
          call ga_set_chunk(g_pr,chunk)
          status = ga_allocate(g_pr)
          call ga_zero(g_pr)
C  Global sync array init
          ndim = 1
          chunk(1) = 1 
          g_stop = ga_create_handle()
          call ga_set_data(g_stop,ndim,dims(2),MT_F_INT)
          call ga_set_chunk(g_stop,chunk(1))
          status =  ga_allocate(g_stop)
          call ga_zero(g_stop)
C  End Global sync array init
          call dntmc_sync(.false.)
       endif
          call ga_pgroup_set_default(idmedium)
          call dntmc_sync(.false.)

C       call ga_pgroup_set_default(idbig)
C       call ga_sync() ! make sure that everyone is at the same place
C       call dntmc_sync(.false.)
C       call dntmc_sync(.false.)
c     endif
C********************************************************
c
c For now I am letting all of the processors run through all of this
c it is probably better in the long term to only have the 0 node
c on the subgroups doing this.

C      iam = myid
C  Change for subgroup only identification
       call dntmc_write_prcfil(idmedium, groupnums, prcfil,ngroups)
       call dntmc_build_prcfil(prcfil2, myid)
       write(6,*)'groupnums for node ',myid,' is ', groupnums(1:ngroups)
       call flush(6) 
C  End Change for subgroup only identification
C*************************************************************
C*************************************************************
C    output files
C*************************************************************

C general output

      write(fprefix,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix))
     &     call errquit('rtdb get file_prefix failed',0,RTDB_ERR)
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(5(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.fort12.',
     &      prcfil(1:4), '.', prcfil2(1:4)

      OPEN(UNIT = 12,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')

C  Added files for configurations
C  Unit 10 ALL
C  Unit 11 Accepted Only
      write(fprefix,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix))
     &     call errquit('rtdb get file_prefix failed',0,RTDB_ERR)
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(5(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.fort10.',
     &      prcfil(1:4), '.', prcfil2(1:4)

      OPEN(UNIT = 10,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')

      write(fprefix,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix))
     &     call errquit('rtdb get file_prefix failed',0,RTDB_ERR)
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(5(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.fort11.',
     &      prcfil(1:4), '.', prcfil2(1:4)

      OPEN(UNIT = 11,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')
C  End Added files for configurations
C      OPEN(UNIT = 12,FILE = file_prefix(1:inp_strlen(file_prefix))
C     &      //'.fort12.'//prcfil(1:4),STATUS =
C     +     'UNKNOWN')

c     OPEN(UNIT = 16,FILE ='10w2-25p888L.'//prcfil(1:4)//'.xmol', 
c    +    STATUS = 'UNKNOWN')

C probability distribution function, permanent
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(5(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.fort3.',
     &      prcfil(1:4),'.',prcfil2(1:4)

      OPEN(UNIT = 3,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')

C      OPEN(UNIT = 3,FILE = file_prefix(1:inp_strlen(file_prefix))
C     &      //'.fort3.'//prcfil(1:4),STATUS =
C     +     'UNKNOWN')

C configuration file, used for restart also, permanent
C      OPEN(UNIT = 9, FILE = '10w2-25p888Lcon.'//prcfil(1:4), 
C     +    STATUS = "UNKNOWN")

C other distribution files, permanent
C      OPEN(UNIT = 10, FILE = '10w2-25p888Lnum.'//prcfil(1:4), 
C     +    STATUS = "UNKNOWN")
C
C      OPEN(UNIT = 11, FILE = '10w2-25p888Lnums.'//prcfil(1:4), 
C     +    STATUS = "UNKNOWN")
C
C      OPEN(UNIT = 13, FILE = '10w2-25p888Lnuma.'//prcfil(1:4), 
C     +    STATUS = "UNKNOWN")

C*************************************************************
C*************************************************************

      kb = 1.987905D-03 
      k = 1.293601055D+40
      e = 1.60219D-19
      PI = 3.141592653589793D+00
      esqrd = e*e  
      
C logical flags
      IACCEPTT = 0
      IACCEPTR = 0
      IACCEPTV = 0
      IACCEPTPR = 0
      IACCEPTPI = 0


C     SPC(NMOL) = PARTICULAR SPECIES -> 
C      ->  (H2SO4=1,H2O=2,H3O+=3,HSO4-=4,H+=5,SO4 -2=6)

C     NATOM(SPC(NMOL)) = # OF ATOMS IN PARTICULAR SPECIES
c needs to be replaced with number of atoms in each species
c as determined in the geometries - generic

      NATOM(1) = 7    !sulfuric acid -> H2SO4
      NATOM(2) = 4    !Dang-Chang water (4sites) -> H2O
      NATOM(3) = 4    !hydronium -> H3O+
      NATOM(4) = 6    !bisulfate -> HSO4-
      NATOM(5) = 1    !proton -> H+
      NATOM(6) = 5    !sulfate -> SO4 -2

C     THESE ARE THE CHARGE, SIGMA, AND EPSILON PARAMETERS
C     FOR THE DANG-CHANG POLARIZABLE WATER MODEL
c this should be sent to an initialization for the potential

      QSE(2,1,1) = 0.0000
      QSE(2,1,2) = 3.2340
      QSE(2,1,3) = 0.1825
      QSE(2,2,1) = 0.5190
      QSE(2,2,2) = 0.0000
      QSE(2,2,3) = 0.0000
      QSE(2,3,1) = 0.5190
      QSE(2,3,2) = 0.0000
      QSE(2,3,3) = 0.0000
      QSE(2,4,1) = -1.0380
      QSE(2,4,2) = 0.0000
      QSE(2,4,3) = 0.0000

C     MASSES FOR WATER
c need to pull this out of the geom in dntmc_getspecinfo
c     M(2,1) = 15.995
c     M(2,2) = 1.008
c     M(2,3) = 1.008
c     M(2,4) = 0.000


C     THIS IS THE INITIAL INPUT
c  Defaults to writing restart information
      CONT1 = 1  !restart: 0 = no restart files, 1 = write restart files
c  Get norestart information from rtdb
      if (.not. rtdb_get(rtdb, 'dntmc:norestart', MT_INT, 1, CONT1))
     &   CONT1 = 1
c
c
c  Defaults to not a restart run
      CONT2 = 0  !is this a restart run?: 0 = no, 1 = yes
c  Get restart information from rtdb
      call util_get_rtdb_state(rtdb, os, oc, or)
      if (or .eq. .true.) then
         CONT2=1
      else
         CONT2=0
      endif
c
c Try to get the geometry related information from the rtdb
c
      call dntmc_getspecinfo(rtdb, nmols, natom, rinit, uspc, nspc,
     &      m)
c
      if (.not.rtdb_get(rtdb, 'dntmc:mcsteps', MT_INT, 1, MCSTPS))
     &    call errquit('dntmc_gibbs: rtdb error for mcsteps', 0, 
     &    RTDB_ERR)


      IF(CONT2.EQ.0)THEN

         LMIN = 1
         LMAX = MCSTPS

      ELSE

C         LMIN = MCSTPS + 1
C         LMAX = 11000000     !assuming previous run is the same length
C  Change for restart handleing 
C  Get info from rtdb
         if (.not. rtdb_get(rtdb, 'dntmc:rsmcstps:'//prcfil(1:4),
     &   MT_INT, 1, rsmcstps))
     &   rsmcstps = 0

         LMIN = rsmcstps
         LMAX = rsmcstps + MCSTPS
C  End Change
      ENDIF
        
C      MPRNT = MCSTPS/10
C  Change for restart handleing
      if(.not. rtdb_get(rtdb, 'dntmc:mprnt', MT_INT, 1, mprnt))
     &  mprnt = 10
C  End Change
C  Define DISP and ALPHA through input directives
C      DISP = 0.03D+00
C      ALPHA = 0.08D+00
      if (.not. rtdb_get(rtdb, 'dntmc:disp', MT_DBL, 1, disp))
     &   disp = 0.03D+00
      if (.not. rtdb_get(rtdb, 'dntmc:alpha', MT_DBL, 1, alpha))
     &   alpha = 0.08D+00
C  End definition of dist and alpha
      ROTF = 0
      AB = 1
      IDUM = -1
C Define a logical for use of a rweight distribution
      if (.not. rtdb_get(rtdb, 'dntmc:rweight', MT_LOG, 1, rweight))
     &   rweight = .false.

C End varible definition

c need to set this based on the nwchem input in dntmc_getspecinfo
c start at 1 for the species numbers
c     USPC(1) = 2
C
C Get from rtdb
c
c      TEMP = 243.D+00
      if (.not.rtdb_get(rtdb, 'dntmc:temp', MT_DBL, 1, TEMP))
     &    call errquit('dntmc_gibbs: rtdb error for temp', 0, RTDB_ERR)


C    INITIALIZE RMIN and RMAX

c
C Get from rtdb
c
c      RMIN = 1.0D+00
c      RMAX = 10.0D+00
      if (.not.rtdb_get(rtdb, 'dntmc:rmin', MT_DBL, 1, RMIN))
     &    call errquit('dntmc_gibbs: rtdb error for rmin', 0, RTDB_ERR)
      if (.not.rtdb_get(rtdb, 'dntmc:rmax', MT_DBL, 1, RMAX))
     &    call errquit('dntmc_gibbs: rtdb error for rmax', 0, RTDB_ERR)
    
      RSIM = RMAX


      DO I = 1, NMOLS

         SPC(I) = USPC(1)

      ENDDO


      WRITE(12,*)NSPC,NMOLS,MCSTPS,MPRNT,TEMP,DISP,ALPHA,ROTF
      WRITE(12,*)AB
      WRITE(12,*)'CONT1 =',CONT1
      WRITE(12,*)'CONT2 =',CONT2     


      WRITE(6,*)NSPC,NMOLS,MCSTPS,MPRNT,TEMP,DISP,ALPHA,ROTF
      WRITE(6,*)AB
      WRITE(6,*)'CONT1 =',CONT1
      WRITE(6,*)'CONT2 =',CONT2


c     WRITE(6,*)'why doesnt this code work'
c     call ga_terminate()
c     call pend
c     STOP

 

C calculate total mass of the system
      DO I = 1, NSPC
         DO J = 1, NATOM(USPC(I))
            MTOT(USPC(I)) = MTOT(USPC(I)) + M(USPC(I),J)
            WRITE(12,*)'M(',USPC(I),',',J,') =',M(USPC(I),J)
         END DO
      END DO
c need to fix this write to be independent of the species number
      WRITE(12,*)'MTOT(2) =',MTOT(2)

C******************************************************************

      IF(CONT2.EQ.0)THEN   !CONT2 = 0 initial run, = 1 restart run

      SEED(1) = 962525439.0

C******************************************************************
C     CREATE INITIAL CONFIGURATION OF WATER MOLECULES
C     CORRESPONDING TO WHICH PARTICULAR NODE THE CODE
C     IS RUNNING ON

c This needs to be on the group ID so that all members of a group
c are doing the same thing.

      NODE = idmedium

      CALL dntmc_CONFIG(NMOLS,NODE,RMAX,rinit)

C  Distribute starting configurations between rmin and rmax
C  Center Cluster
C center of mass
      CALL dntmc_CLUSTERCMASS(NMOLS)

C recenter
      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO
      
C find outermost molecule and vector
      DO I = 1, NMOLS
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DISTcalc(I) = RX*RX + RY*RY + RZ*RZ
      END DO
      DIST = DISTcalc(1)
      molmax = 1 
      DO I = 2, NMOLS
         IF (DISTcalc(I) .gt. DIST) then
            DIST = DISTcalc(I)
            molmax = I
         endif
      enddo
      NODE = -1
      do i=1, ngroups
        if (groupnums(i) .eq. idmedium) then
          NODE = i
        endif
      enddo
      DIST = DSQRT(DIST)
      RX = RCM(molmax,1)/DIST
      RY = RCM(molmax,2)/DIST
      RZ = RCM(molmax,3)/DIST

C  Manual Scanning
      if(.not. rtdb_get(rtdb,'dntmc:node',MT_INT,1,node)) then
      NODE = 1 
      endif
C
      DIST = RMIN - DIST
C Manual Scanning
C      DIST = DBLE(NODE)*(RMAX-RMIN)/DBLE(ngroups+1) + DIST
      DIST = DBLE(NODE)*(RMAX-RMIN)/100.00D+00 + DIST
C
      if (NMOLS .eq. 1) then
        DIST = 0.00D+00
      else
        DIST = (DBLE(NMOLS)/DBLE(NMOLS -1)) * DIST
      endif
      NODE = idmedium
C Translate by constant
      DO I = 1, NATOM(SPC(molmax))
            R(molmax,I,1) = R(molmax,I,1)+(RX * DIST)
            R(molmax,I,2) = R(molmax,I,2)+(RY * DIST)
            R(molmax,I,3) = R(molmax,I,3)+(RZ * DIST)
      END DO

C Recenter
      CALL dntmc_CLUSTERCMASS(NMOLS)
      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

C set rsim if rweight method
      if(rweight) then
      DO I = 1, NMOLS
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DISTcalc(I) = RX*RX + RY*RY + RZ*RZ
      END DO
      DIST = DISTcalc(1)
      DO I = 2, NMOLS
         IF (DISTcalc(I) .gt. DIST) then
            DIST = DISTcalc(I)
         endif
      enddo
      rsim = DSQRT(DIST)
      endif
C end of rsim setting

C End Distribute configurations 

      WRITE(6,*)'made it through the initial config setup'
      WRITE(6,*)'on group =',prcfil(1:4)

      WRITE(12,*)'made it through the initial config setup'
      WRITE(12,*)'on group =',prcfil(1:4)


C*********************************************************

      ELSE       !CONT2 = 1 restart run
C  Change for restart from rtdb
        write(6,*)'Reading Restart Information for group ',prcfil(1:4)
        call flush(6)
        write(12,*)'Reading Restart Information for group ',prcfil(1:4)
        call flush(12)
        call dntmc_readrestart(rtdb, uu, uavr, usqrd,
     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +  iacceptpi, rsh, prcfil, nmols, rsmcstps)
        write(6,*)'Done with restart setup for group ',prcfil(1:4)
        call flush(6)
        write(12,*)'Done with restart setup for group ',prcfil(1:4)
        call flush(12)
C  End Change for restart

      ENDIF   !endif for cont2 = 1 running from restart

C******************************************************
C     THESE INDICES CONTROL WHICH MOLECULES ARE PROBES
c need to figure this out with Shawn
      DO I = 1, NMOLS - 1
         LNDX(I) = 0     
      END DO

      LNDX(NMOLS) = 1

C*********************************************************
C*********************************************************

      LAMBDA = 1.D-02
      
      ABNDX(1) = -1.D+00      ! B ENSEMBLE
      ABNDX(2) = 1.D+00      ! A ENSEMBLE

C     THESE INDICES ARE FOR THE (H+ SO4 2=)(H+ H+) PROBE INTERACTION = 1
c guessing that this needs to be made generic instead of specific
c need to figure out what they mean

      LNDXX(6,5) = 1      !SO4-2  H+
      LNDXX(5,6) = 1      !H+ SO4-2
      LNDXX(6,2) = 0      !SO4-2 H2O
      LNDXX(2,6) = 0      !H2O SO4-2
      LNDXX(5,2) = 0      !H+ H2O
      LNDXX(2,5) = 0      !H2O H+
      LNDXX(5,5) = 1      !H+ H+
      LNDXX(6,6) = 0      !SO4-2 SO4-2
      LNDXX(2,2) = 0      !H2O H2O
      

C *******************************************************************

c are these dependent on the total number of molecules or the
c total for each species?
      NMIN = NMOLS - 4

      NCURR = NMOLS

      NMAX = NMOLS

c     initialization information for rsim and rmax    

      RSIMSQ = RSIM*RSIM

      RMAXSQ = RMAX*RMAX

      WRITE(12,*)'RMAXSQ =',RMAXSQ
      WRITE(12,*)'RSIMSQ =',RSIMSQ
      WRITE(12,*)'TEMPERATURE =',TEMP,' K'
      WRITE(12,*)'RMIN =',RMIN,' ANGSTROMS'
      WRITE(12,*)'RMAX =',RMAX,' ANGSTROMS'
      WRITE(12,*)'The initial simulation radius =',RSIM,' ANGSTROMS' 

      WRITE(12,*)'NCURR =',NCURR

C Skip check for restart
      IF(CONT2 .eq. 1) goto 101
C Skip check for rweight
      IF(rweight) goto 101
C End Skip


      DO I = 1, NCURR
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DIST = RX*RX + RY*RY + RZ*RZ

         IF(DIST.GE.RSIMSQ) goto 99
      END DO

C  Original Check
C      DO I = 1, NCURR
C         DO J = 1, NATOM(SPC(I))
C              RX = R(I,J,1)
C              RY = R(I,J,2)
C              RZ = R(I,J,3)
C
C              DIST = RX*RX + RY*RY + RZ*RZ
C
C              IF(DIST.GE.RSIMSQ) goto 99
C
C         END DO
C      END DO

      GOTO 101

  99  CONTINUE
      WRITE(12,*)'dist is ',DSQRT(DIST),' mol ',i,' atom ',j
C should never get here unless the algorithm is messed up
      WRITE(12,*)'molecule',I,' is out of the simulation volume'
      call errquit('dntmc_gibbs: should not get to 99', 1, UNKNOWN_ERR)

 101  CONTINUE
      
C grid of pairs of interactions squared
c is this for the Liam-Dang potential?
      DO I = 1, NSPC
      
      DO II = 1, NATOM(USPC(I))
      
      DO JJ = 1, NATOM(USPC(I))
      
      QSQRD(USPC(I),USPC(I),II,JJ) = QSE(USPC(I),II,1)*QSE(USPC(I),JJ,1)
      QSQRD(USPC(I),USPC(I),JJ,II) = QSQRD(USPC(I),USPC(I),II,JJ)
               
      SGM = 0.5D+00*(QSE(USPC(I),II,2) + QSE(USPC(I),JJ,2))
      
      SGMSQ(USPC(I),USPC(I),II,JJ) = SGM*SGM
      SGMSQ(USPC(I),USPC(I),JJ,II) = SGMSQ(USPC(I),USPC(I),II,JJ)
      
      EPSLN(USPC(I),USPC(I),II,JJ) = DSQRT(QSE(USPC(I),II,3)
     +*QSE(USPC(I),JJ,3))
      EPSLN(USPC(I),USPC(I),JJ,II) = EPSLN(USPC(I),USPC(I),II,JJ)      

      END DO
      
      END DO
      
      END DO      
      

      DO I = 1, NSPC - 1
      
      DO J = I + 1, NSPC
         
      DO II = 1, NATOM(USPC(I))
            
      DO JJ = 1, NATOM(USPC(J))
      
      QSQRD(USPC(I),USPC(J),II,JJ) = QSE(USPC(I),II,1)*QSE(USPC(J),JJ,1)
      QSQRD(USPC(J),USPC(I),JJ,II) = QSQRD(USPC(I),USPC(J),II,JJ)
      
      SGM = 0.5D+00*(QSE(USPC(I),II,2) + QSE(USPC(J),JJ,2))
               
      SGMSQ(USPC(I),USPC(J),II,JJ) = SGM*SGM
      SGMSQ(USPC(J),USPC(I),JJ,II) = SGMSQ(USPC(I),USPC(J),II,JJ)

      EPSLN(USPC(I),USPC(J),II,JJ) = DSQRT(QSE(USPC(I),II,3)
     +*QSE(USPC(J),JJ,3))
     
      EPSLN(USPC(J),USPC(I),JJ,II) = EPSLN(USPC(I),USPC(J),II,JJ)

      END DO
               
      END DO
            
      END DO
         
      END DO
            
C******************************************************************

C center of mass
      CALL dntmc_CLUSTERCMASS(NCURR)

C recenter
      DO I = 1, NCURR
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

C out of the sphere after recentering
C Skip check for restart
      IF(CONT2 .eq. 1) goto 456
C skip chesk if rweight
      if(rweight) goto 455
C End Skip

      DO I = 1, NCURR
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DIST = RX*RX + RY*RY + RZ*RZ

         IF(DIST.GT.RSIMSQ)THEN

         WRITE(12,*)'DIST =',DIST

         WRITE(12,*)'MOLECULE ',I,' IS OUT OF SIMULATION VOLUME!'
         call errquit('dntmc_gibbs: should not get here',
     +                 1, UNKNOWN_ERR)

         ENDIF

      END DO

C  Original Code
C      DO I = 1, NCURR
C
C         RX = R(I,1,1)
C         RY = R(I,1,2)
C         RZ = R(I,1,3)
C
C         DIST = RX*RX + RY*RY + RZ*RZ
C
C         IF(DIST.GT.RSIMSQ)THEN
C
C         WRITE(12,*)'DIST =',DIST
C
C         WRITE(12,*)'MOLECULE ',I,' IS OUT OF SIMULATION VOLUME!'
C         call errquit('dntmc_gibbs: should not get here', 
C     +                 1, UNKNOWN_ERR)
C
C         ENDIF
C
C      END DO
C Skip to if rweight 
 455  CONTINUE

c234567      
      WRITE(12,*)'The initial cluster coordinates for a fixed center'
      WRITE(12,*)'of mass at the origin are:'

      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            WRITE(12,*)(R(I,J,IR),IR=1,3)
            RI(I,J,1) = R(I,J,1)
            RI(I,J,2) = R(I,J,2)
            RI(I,J,3) = R(I,J,3)
         END DO
      END DO
      
C check to make sure center of mass is origin
      CALL dntmc_CLUSTERCMASS(NCURR)


C calculate the energy of the configuration
c      CALL dntmc_PSUM(NCURR,UU)
C  Added sync for groups
C           call dntmc_sync(.false.)
C  End added sync
      call dntmc_setcoords(rtdb, r, nmols, natom, spc)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          write(12,*) 'Energy did not converge'
          uu = -0.5                        ! set the energy arbitrarily high - H atom
        else
          write(12,*) 'Energy converged on the second try'
          if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, uu))
     $      call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
        endif
      else 
        if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, uu))
     $     call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
      endif
c      if (.not. task_energy(rtdb)) call errquit
c     $     ('dntmc_gibbs: energy failed',0, GEOM_ERR)
c      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, uu))
c     $     call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
c
c convert to kcal/mol
c
      uu = uu*627.5095D+00
C Skip to if restart
 456  CONTINUE
      WRITE(12,*)'INITIAL UU = ',UU,' kcal/mole'

C calculate number of atoms
      DO I = 1, NMOLS
         NATOT = NATOT + NATOM(SPC(I))
c	 PRINT*,'NATOT =',NATOT	 
      END DO
      
      INAME = 97 - 1
      
      WRITE(12,*)'LMIN =',LMIN,'  LMAX =',LMAX
       
      WRITE(6,*)'made it through initialization'

      IICOUNT = 0

C beginning of MC loop
      DO L = LMIN , LMAX

         ACCEPTT = .FALSE. ! translate
         ACCEPTR = .FALSE. ! rotation
         ACCEPTV = .FALSE. ! volume
         ACCEPTPR = .FALSE. !????
         ACCEPTPI = .FALSE.
         
C do one move, calculate energy and return the acceptance ratios and energy
      if(rweight) then
         CALL dntmc_MOVE_RWEIGHT(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
     +             ,ACCEPTPI,UU, rtdb, nmols)
      else
         CALL dntmc_MOVE(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
     +             ,ACCEPTPI,UU, rtdb, nmols)
      endif
         CALL dntmc_ROOTMEANSQUARE(L,NCURR,NSPC)

         IF(ACCEPTT) IACCEPTT = IACCEPTT + 1
         IF(ACCEPTR) IACCEPTR = IACCEPTR + 1
         IF(ACCEPTV) IACCEPTV = IACCEPTV + 1
         IF(ACCEPTPR) IACCEPTPR = IACCEPTPR + 1
         IF(ACCEPTPI) IACCEPTPI = IACCEPTPI + 1

         ARATIOTT = (DBLE(IACCEPTT)/DBLE(L))*100.D+00
         ARATIOTR = (DBLE(IACCEPTR)/DBLE(L))*100.D+00
         ARATIOTV = (DBLE(IACCEPTV)/DBLE(L))*100.D+00
         ARATIOTPR = (DBLE(IACCEPTPR)/DBLE(L))*100.D+00
         ARATIOTPI = (DBLE(IACCEPTPI)/DBLE(L))*100.D+00

C calculate average of stuff
         CALL dntmc_AVERAGE(L,TEMP,UU,UAVR,Cv,RSIMAVR)
c calculate distribution of different radii of the sphere
         CALL dntmc_NUMBERHIST(L)

       if(rweight) then
       RINC = (RMAX - RMIN)/100.D+00
       RSH = RMIN
       DO I = 1, NOB
          RSH = RSH + RINC
          RB(I) = RSH
       ENDDO
       endif

         CALL dntmc_NUMBERAVERAGE(L,temp,rweight)
C***************************************************************

      IF (MOD(L,MPRNT).EQ.0)THEN

      WRITE(12,*)'*******************************'
      WRITE(12,*)'THIS IS THE END OF STEP = ',L
      WRITE(12,*)'the number of molecules currently'
      WRITE(12,*)'in the simulation are',NCURR
      WRITE(12,*)'ARATIOT = ',ARATIOTT,' %'
      WRITE(12,*)'ARATIOR = ',ARATIOTR,' %'
      WRITE(12,*)'ARATIOV = ',ARATIOTV,' %'
      WRITE(12,*)'ARATIOPR = ',ARATIOTPR,' %'
      WRITE(12,*)'ARATIOPI = ',ARATIOTPI,' %'
      WRITE(12,*)'UNEW =',UU,' kcal/mole ' 
      WRITE(12,*)'RSIM =',RSIM,' ANGSTROMS'
      WRITE(12,*)'<RSIM> =',RSIMAVR,' ANGSTROMS'
      WRITE(12,*)'<U> =',UAVR,' kcal/mole' 
      WRITE(12,*)'<Cv> =',Cv,' kcal/mole-K'
      WRITE(12,*)'*******************************'      


C*******************************************************************
C*******************************************************************
C     writing occasional restart files
      lnow = L + 1
      IF(CONT1.EQ.1)THEN
      if (ga_pgroup_nodeid(idmedium) .eq. 0) then
        call dntmc_writefiles(uu, uavr, usqrd,
     +    rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +    iacceptpi, rsh, prcfil, nmols, lnow, fprefix)
        write(6,*)' Restart files written for group ',prcfil(1:4)
        call flush(6)
        write(12,*)' Restart files written for group ',prcfil(1:4)
        call flush(12)
      endif
C  End Restart changes          
      ENDIF     !this endif is for writing restart files

C*******************************************************************
C*******************************************************************
C     THIS SECTION IS FOR CHECKING DISTRIBUTIONS FROM ALL
C     NODES AND ENSURES THAT THE MAXIMUM VARIANCE IS SMALL 
C     (<15%)


C     STORE RADIAL GRID RSH INTO RB

       RINC = (RMAX - RMIN)/100.D+00
       RSH = RMIN
       DO I = 1, NOB
          RSH = RSH + RINC
          RB(I) = RSH
       ENDDO

C**********************************************************************
C    THIS SECTION STORES THE LOCAL NUMBERAVR INTO THE GLOBAL ARRAY
C      call ga_pgroup_set_default(idzero)
      if (ga_pgroup_nodeid(idmedium).eq.0) then
         call ga_pgroup_set_default(idzero)
C  Global sync mutex
         ld(1) = 1
         ld(2) = 1
         call ga_put(g_stop,ga_nodeid()+1,ga_nodeid()+1,1,1,ld(1),ld)
         call ga_lock(0)
         write(6,*)' I made it into mutex', ga_nodeid()
         call flush(6)
         trigger = .true.
         do while (trigger)
         ld(1) = ngroups
         ld(2) = 1
         call ga_get(g_stop,1,ngroups,1,1,stoplist,ld)
         ld(1) = 0
         do i=1, ngroups
           ld(1) = ld(1) + stoplist(i)
         end do
         if (ld(1) .lt. ngroups) then
            trigger = .true.
         else
            trigger = .false.
         endif
         end do
         call ga_unlock(0)
         write(6,*)' I made it out of mutex', ga_nodeid()
         call flush(6)
C  Global sync mutex end
         ilo = 1
         ihi = NOB
         jlo = ga_nodeid() + 1
         jhi = ga_nodeid() + 1

c      write(6,*)'jlo =',jlo

         ld(1)  = 100   !must be physical dimension of local array
         ld(2)  = 1
         call ga_put(g_pr,ilo,ihi,jlo,jhi,NUMBERAVR,ld)
C         call ga_sync
          call dntmc_sync(.false.)
C         call ga_pgroup_sync(idbig)
         jlo = 1
         jhi = ga_nnodes()

c      write(6,*)'jhi =',jhi
c
c everyone gets the full results

         call ga_get(g_pr,ilo,ihi,jlo,jhi,TESTME,ld)
         call dntmc_sync(.false.)
C*********************************************************************
C     this section is for testing...

         RSH = RMIN
         DO I = 1, NOB
            RSH = RSH + RINC
C            WRITE(6,*)RSH,TESTME(I,ga_nodeid()+1)
         ENDDO

C**********************************************************************
C    this section is averages the distributions from all nodes
C    and finds the maximum variance and tests it against
C    a predefined tolerance in the logical below...
C  Added nda_old array
         CALL dntmc_AVERAGEDIST(NOB,ga_nnodes(),PERCENT,nda_old)
         write(6,*) 'ga_nnodes and percent in idzero ',
     &     ga_nnodes(), percent
C  Global sync reinit
      call ga_zero(g_stop)
C  End global sync reinit
      endif  ! idzero section
      call ga_pgroup_set_default(idmedium)
C      call dntmc_sync(.false.)
C      call ga_sync()
C           call dntmc_sync(.false.)
C      call ga_brdcst(msg_dntmc0+MSGINT, PERCENT, mitob(1), 0)
      call ga_brdcst(msg_dntmc0+MSGDBL, PERCENT, mitob(1), 0)
      call ga_brdcst(msg_dntmc1+MSGDBL, NDA, mdtob(100), 0)
      call ga_brdcst(msg_dntmc2+MSGDBL, NDSIG, mdtob(100), 0)
      call dntmc_sync(.false.)
C  Added save of old NDA array
         Do i=1, NOB
         nda_old(i) = NDA(i)
         enddo
C  End added old save

      write(6,*) 'percent is ', percent, ga_nodeid(), myid
C      IF(PERCENT.LT.30.D+00)THEN
       IF(PERCENT.LT.(100.D+00/DBLE(NOB)))THEN

          WRITE(6,*)'the average distribution has converged'
          WRITE(12,*)'the average distribution has converged'
          WRITE(6,*)'PERCENT =',PERCENT
          WRITE(12,*)'PERCENT =',PERCENT

          GOTO 150

      ELSE

          WRITE(6,*)'the average distribution has not converged yet'
          WRITE(12,*)'the average distribution has not converged yet'
          WRITE(6,*)'PERCENT =',PERCENT
          WRITE(12,*)'PERCENT =',PERCENT

c gather best estimate of R
          CALL dntmc_RANGER(NOB,Rest)

          WRITE(6,*)'best estimate of Rcut =',Rest
          WRITE(12,*)'best estimate of Rcut =',Rest      

      ENDIF


C*******************************************************************
C*******************************************************************
C    this endif is for periodic checks for small variance on
C    the distributions, printouts, etc.    

      ENDIF     

c**************************************************************
C**************************************************************
C     THIS ENDDO IS FOR THE MAIN MC LOOP     

      ENDDO

c**************************************************************
C**************************************************************

      WRITE(12,*)'did not converge average distribution'
      WRITE(6,*)'did not converge average distribution'

C    WRITING FINAL AVERAGE DISTRIBUTION WITH ERROR BARS
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(5(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.fort30.',
     &      prcfil(1:4),'.',prcfil2(1:4)

      OPEN(UNIT = 30,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS ='REPLACE')

C      OPEN(UNIT = 30,FILE = file_prefix(1:inp_strlen(file_prefix))
C     &      //'.fort30.'//prcfil(1:4),STATUS =
C     +     'UNKNOWN')

C      DO I = 1, NOB - 1
       DO I = 1, NOB
         WRITE(30,*)RB(I),NDA(I),NDSIG(I)
      END DO

      CLOSE(30)

      GOTO 180
     
c**************************************************************
C**************************************************************
c**************************************************************
c**************************************************************
C    THIS IS THE BEGINNING OF A NEW MONTE CARLO RUN RESULTING 
C    FROM A BETTER ESTIMATE OF Rcut*

 150  CONTINUE   

      WRITE(6,*)'going into dntmc_RANGER again for checking'
      WRITE(12,*)'going into dntmc_RANGER again for checking'            

      CALL dntmc_RANGER(NOB,Rest)

      WRITE(6,*)'best estimate of Rcut =',Rest
      WRITE(12,*)'best estimate of Rcut =',Rest      
      
      IF(Rest.LT.DBLE(NMOLS))THEN
         Rest = DBLE(NMOLS) + 2.D+00
      ENDIF

      RMAX = Rest + 3.D+00   !3 angstroms beyond the best estimate

      WRITE(6,*)'new value of RMAX =',RMAX
      WRITE(12,*)'new value of RMAX =',RMAX 

C     STORE NEW RADIAL GRID RSH INTO RB

       RINC = (RMAX - RMIN)/100.D+00
       RSH = RMIN
       DO I = 1, NOB
          RSH = RSH + RINC
          RB(I) = RSH
       ENDDO

C    re-initialize global array to zero
C      call ga_pgroup_set_default(idzero)
      if(ga_pgroup_nodeid(idmedium).eq.0) then
         call ga_pgroup_set_default(idzero)
         CALL ga_zero(g_pr)
         call dntmc_sync(.false.)
      endif
C      call dntmc_sync(.false.)
      call ga_pgroup_set_default(idmedium)
      call dntmc_sync(.false.)
C      call ga_sync()
C           call dntmc_sync(.false.)

C******************************************************************
C     CREATE NEW INITIAL CONFIGURATION OF WATER MOLECULES
C     CORRESPONDING TO WHICH PARTICULAR NODE THE CODE
C     IS RUNNING ON
C     THIS IS THE SECOND MONTE CARLO LOOP USING A BETTER Rcut*

      NODE = idmedium

      CALL dntmc_CONFIG(NMOLS,NODE,RMAX,rinit)

C*********************************************************
c     initialization information for rsim and rmax   

      RSIM = RMAX 

      RSIMSQ = RSIM*RSIM

      RMAXSQ = RMAX*RMAX

      WRITE(6,*)'made it through new configuration initialization'
      WRITE(12,*)'made it through new configuration initialization'

      WRITE(12,*)'RMAXSQ =',RMAXSQ
      WRITE(12,*)'RSIMSQ =',RSIMSQ
      WRITE(12,*)'TEMPERATURE =',TEMP,' K'
      WRITE(12,*)'RMIN =',RMIN,' ANGSTROMS'
      WRITE(12,*)'RMAX =',RMAX,' ANGSTROMS'
      WRITE(12,*)'The initial simulation radius =',RSIM,' ANGSTROMS' 
      WRITE(12,*)'NCURR =',NCURR
      WRITE(6,*)'RMAXSQ =',RMAXSQ
      WRITE(6,*)'RSIMSQ =',RSIMSQ
      WRITE(6,*)'TEMPERATURE =',TEMP,' K'
      WRITE(6,*)'RMIN =',RMIN,' ANGSTROMS'
      WRITE(6,*)'RMAX =',RMAX,' ANGSTROMS'
      WRITE(6,*)'The initial simulation radius =',RSIM,' ANGSTROMS' 
      WRITE(6,*)'NCURR =',NCURR

C    Checking to see if initial configuration is in 
C    new simulation volume

      DO I = 1, NCURR
         DO J = 1, NATOM(SPC(I))
              RX = R(I,J,1)
              RY = R(I,J,2)
              RZ = R(I,J,3)
              DIST = RX*RX + RY*RY + RZ*RZ
              IF(DIST.GE.RSIMSQ) goto 160
         END DO
      END DO

      GOTO 170

 160  CONTINUE

C you should never get here
      WRITE(12,*)'molecule',I,' is out of the simulation volume'
      call errquit('dntmc_gibbs: should not get to 160', 1, UNKNOWN_ERR)

 170  CONTINUE

      CALL dntmc_CLUSTERCMASS(NCURR)

      DO I = 1, NCURR
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

      DO I = 1, NCURR

         RX = R(I,1,1)
         RY = R(I,1,2)
         RZ = R(I,1,3)

         DIST = RX*RX + RY*RY + RZ*RZ

         IF(DIST.GT.RSIMSQ)THEN
            WRITE(12,*)'DIST =',DIST
            WRITE(12,*)'MOLECULE ',I,' IS OUT OF SIMULATION VOLUME!'
            call errquit('dntmc_gibbs: should not get to 160', 
     +                    1, UNKNOWN_ERR)
         ENDIF

      END DO
c234567      
      WRITE(12,*)'The initial cluster coordinates for a fixed center'
      WRITE(12,*)'of mass at the origin are:'

      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            WRITE(12,*)(R(I,J,IR),IR=1,3)
            RI(I,J,1) = R(I,J,1)
            RI(I,J,2) = R(I,J,2)
            RI(I,J,3) = R(I,J,3)
         END DO
      END DO
      
      CALL dntmc_CLUSTERCMASS(NCURR)

c      CALL dntmc_PSUM(NCURR,UU)
      call dntmc_setcoords(rtdb, r, nmols, natom, spc)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          write(12,*) 'Energy did not converge'
          uu = -0.5                        ! set the energy arbitrarily high - H atom
        else
          write(12,*) 'Energy converged on the second try'
          if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, uu))
     $      call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
        endif
      else 
        if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, uu))
     $     call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
      endif
c      if (.not. task_energy(rtdb)) call errquit
c     $     ('dntmc_gibbs: energy failed',0, GEOM_ERR)
c      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, uu))
c     $     call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
     
c
c convert to kcal/mol
c
      uu = uu*627.5095D+00
      WRITE(12,*)'2nd MC LOOP INITIAL UU = ',UU,' kcal/mole'
      WRITE(6,*)'2nd MC LOOP INITIAL UU = ',UU,' kcal/mole'

C    RE-INITIALIZE ACCEPTANCE RATIOS

      IACCEPTT = 0
      IACCEPTR = 0
      IACCEPTV = 0
      IACCEPTPR = 0
      IACCEPTPI = 0

      DO L = LMIN , LMAX

         ACCEPTT = .FALSE.
         ACCEPTR = .FALSE.
         ACCEPTV = .FALSE.
         ACCEPTPR = .FALSE.
         ACCEPTPI = .FALSE.

      if(rweight) then
         CALL dntmc_MOVE_RWEIGHT(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
     +             ,ACCEPTPI,UU, rtdb, nmols)
      else
         CALL dntmc_MOVE(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
     +             ,ACCEPTPI,UU, rtdb, nmols)
      endif

         CALL dntmc_ROOTMEANSQUARE(L,NCURR,NSPC)

         IF(ACCEPTT) IACCEPTT = IACCEPTT + 1
         IF(ACCEPTR) IACCEPTR = IACCEPTR + 1
         IF(ACCEPTV) IACCEPTV = IACCEPTV + 1
         IF(ACCEPTPR) IACCEPTPR = IACCEPTPR + 1
         IF(ACCEPTPI) IACCEPTPI = IACCEPTPI + 1

         ARATIOTT = (DBLE(IACCEPTT)/DBLE(L))*100.D+00
         ARATIOTR = (DBLE(IACCEPTR)/DBLE(L))*100.D+00
         ARATIOTV = (DBLE(IACCEPTV)/DBLE(L))*100.D+00
         ARATIOTPR = (DBLE(IACCEPTPR)/DBLE(L))*100.D+00
         ARATIOTPI = (DBLE(IACCEPTPI)/DBLE(L))*100.D+00

         CALL dntmc_AVERAGE(L,TEMP,UU,UAVR,Cv,RSIMAVR)
         CALL dntmc_NUMBERHIST(L)

       if(rweight) then
       RINC = (RMAX - RMIN)/100.D+00
       RSH = RMIN
       DO I = 1, NOB
          RSH = RSH + RINC
          RB(I) = RSH
       ENDDO
       endif

         CALL dntmc_NUMBERAVERAGE(L,temp,rweight)
C***************************************************************

      IF (MOD(L,MPRNT).EQ.0)THEN

      WRITE(12,*)'*******************************'
      WRITE(12,*)'THIS IS THE END OF STEP = ',L
      WRITE(12,*)'the number of molecules currently'
      WRITE(12,*)'in the simulation are',NCURR
      WRITE(12,*)'ARATIOT = ',ARATIOTT,' %'
      WRITE(12,*)'ARATIOR = ',ARATIOTR,' %'
      WRITE(12,*)'ARATIOV = ',ARATIOTV,' %'
      WRITE(12,*)'ARATIOPR = ',ARATIOTPR,' %'
      WRITE(12,*)'ARATIOPI = ',ARATIOTPI,' %'
      WRITE(12,*)'UNEW =',UU,' kcal/mole ' 
      WRITE(12,*)'RSIM =',RSIM,' ANGSTROMS'
      WRITE(12,*)'<RSIM> =',RSIMAVR,' ANGSTROMS'
      WRITE(12,*)'<U> =',UAVR,' kcal/mole' 
      WRITE(12,*)'<Cv> =',Cv,' kcal/mole-K'
      WRITE(12,*)'*******************************'

C*******************************************************************
C*******************************************************************
C     writing occasional restart files
      lnow = L + 1
      IF(CONT1.EQ.1)THEN
      if (ga_pgroup_nodeid(idmedium) .eq. 0) then
        call dntmc_writefiles(uu, uavr, usqrd,
     +    rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +    iacceptpi, rsh, prcfil, nmols, lnow, fprefix)
        write(6,*)' Restart files written for group ',prcfil(1:4)
        call flush(6)
        write(12,*)' Restart files written for group ',prcfil(1:4)
        call flush(12)
      endif
C  End Restart changes          
      ENDIF     !this endif is for writing restart files

C*******************************************************************
C*******************************************************************

C**********************************************************************
C      call ga_pgroup_set_default(idzero)
      if (ga_pgroup_nodeid(idmedium).eq.0) then
         call ga_pgroup_set_default(idzero)
C  Global sync mutex
         ld(1) = 1
         ld(2) = 1
         call ga_put(g_stop,ga_nodeid()+1,ga_nodeid()+1,1,1,ld(1),ld)
         call ga_lock(0)
         write(6,*)' I made it into mutex', ga_nodeid()
         call flush(6)
         trigger = .true.
         do while (trigger)
         ld(1) = ngroups
         ld(2) = 1
         call ga_get(g_stop,1,ngroups,1,1,stoplist,ld)
         ld(1) = 0
         do i=1, ngroups
           ld(1) = ld(1) + stoplist(i)
         end do
         if (ld(1) .lt. ngroups) then
            trigger = .true.
         else
            trigger = .false.
         endif
         end do
         call ga_unlock(0)
         write(6,*)' I made it out of mutex', ga_nodeid()
         call flush(6)
C  Global sync mutex end
         ilo = 1
         ihi = NOB
         jlo = ga_nodeid() + 1
         jhi = ga_nodeid() + 1

c      write(6,*)'jlo =',jlo

         ld(1)  = 100   !must be physical dimension of local array
         ld(2)  = 1

         call ga_put(g_pr,ilo,ihi,jlo,jhi,NUMBERAVR,ld)
C         call ga_sync()
           call dntmc_sync(.false.)
C          call ga_pgroup_sync(idbig)

         jlo = 1
         jhi = ga_nnodes()

c      write(6,*)'jhi =',jhi

         call ga_get(g_pr,ilo,ihi,jlo,jhi,TESTME,ld)
         call dntmc_sync(.false.)
C*********************************************************************
C  Added nda_old array
         CALL dntmc_AVERAGEDIST(NOB,ga_nnodes(),PERCENT,nda_old)
C  Global sync reinit
      call ga_zero(g_stop)
C  End global sync reinit
      endif  ! idzero section
      call ga_pgroup_set_default(idmedium)
C      call dntmc_sync(.false.)
C      call ga_sync()
C           call dntmc_sync(.false.)
C      call ga_brdcst(msg_dntmc0+MSGINT, PERCENT, mitob(1), 0)
      call ga_brdcst(msg_dntmc4+MSGDBL, PERCENT, mdtob(1), 0)
      call ga_brdcst(msg_dntmc5+MSGDBL, NDA, mdtob(100), 0)
      call ga_brdcst(msg_dntmc6+MSGDBL, NDSIG, mdtob(100), 0)
      call dntmc_sync(.false.)
C  Added save of old NDA array
         Do i=1, NOB
         nda_old(i) = NDA(i)
         enddo
C  End added old save

      WRITE(6,*)'PERCENT =',PERCENT
      WRITE(12,*)'PERCENT =',PERCENT

      ENDIF

c**************************************************************
C**************************************************************
C     THIS ENDDO IS FOR THE NEW MAIN MC LOOP 
       IF(PERCENT.LT.(50.D+00/DBLE(NOB)))THEN

          WRITE(6,*)'the average distribution has converged'
          WRITE(12,*)'the average distribution has converged'
          WRITE(6,*)'PERCENT =',PERCENT
          WRITE(12,*)'PERCENT =',PERCENT

          GOTO 175

      ELSE

          WRITE(6,*)'the average distribution has not converged yet'
          WRITE(12,*)'the average distribution has not converged yet'
          WRITE(6,*)'PERCENT =',PERCENT
          WRITE(12,*)'PERCENT =',PERCENT
      ENDIF

      ENDDO
c**************************************************************
C**************************************************************
C Added end of DO loop destination
 175  CONTINUE
C End added destination

C    WRITING FINAL AVERAGE DISTRIBUTION WITH ERROR BARS
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(5(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.fort24.',
     &      prcfil(1:4),'.',prcfil2(1:4)

      OPEN(UNIT = 24,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')

C      OPEN(UNIT = 24,FILE = file_prefix(1:inp_strlen(file_prefix))
C     &      //'.fort24.'//prcfil(1:4),STATUS =
C     +     'UNKNOWN')

C      DO I = 1, NOB - 1
       DO I =1, NOB

         WRITE(24,*)RB(I),NDA(I),NDSIG(I)

      END DO

      CLOSE(24)

c**************************************************************
C**************************************************************
C    DIDNT CONVERGE AVERAGE DISTRIBUTION
C    BYPASSING 2nd MONTE CARLO LOOP

 180  CONTINUE

c**************************************************************
C**************************************************************

      CALL dntmc_CLUSTERCMASS(NCURR)
c      CALL dntmc_ANGLEAVR(L)

      DO I = 1, 3
         WRITE(12,*)'RCLCM(',I,') =',RCLCM(I)
      END DO

C*******************************************************************
C     WRITE THE HISTOGRAM DATA
C
C     NOB is defind above = number of bins

      RINC = (RMAX - RMIN)/100.D+00
      RSH = RMIN
C      DO I = 1, NOB - 1
      DO I =1, NOB
         RSH = RSH + RINC
         WRITE(3,*)RSH,NUMBERAVR(I)
      END DO

      CLOSE(3)

C*******************************************************************
C     writing restart files
      IF(CONT1.EQ.1)THEN
C  Writing final restart files
C      call ga_pgroup_set_default(idmedium)
      if (ga_nodeid() .eq. 0) then
        call dntmc_writefiles(uu, uavr, usqrd,
     +    rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +    iacceptpi, rsh, prcfil, nmols, L, fprefix)
        write(6,*)' Restart files written for group ',prcfil(1:4)
        call flush(6)
        write(12,*)' Restart files written for group ',prcfil(1:4)
        call flush(12)
      endif
C      call dntmc_sync(.false.)
C  End Restart changes          
          
      ENDIF   !this endif is for writing the final restart files...

      CLOSE(12)
C  Added files for configurations 10 = all 11 = accepted
      CLOSE(10)
      CLOSE(11)
C  End Added files for configurations
      WRITE(6,*)'the end - Ive written everything on each node'

c*************************************
      if (ga_nodeid().eq.0) then
C        call ga_pgroup_set_default(idmedium)
        call ga_pgroup_set_default(idzero)
        status = ga_destroy(g_pr)
C  Global sync array destroy
        status = ga_destroy(g_stop)
C  End Global sync arrey destroy
      else
C        call ga_pgroup_set_default(idmedium)
      endif
      call ga_pgroup_set_default(idmedium)
C  Global sync mutex destroy
      status = ga_destroy_mutexes()
C  End Global sync mutex destroy
      status = ga_destroy(g_proc)
      call util_sgend(rtdb)  ! finalize subgroup stuff
      write(6,*) 'leaving dntmc_gibbs', myid

 9945 Continue ! How to instantly exit code
      dntmc_gibbs = .true.
      write(6,*) 'leaving dntmc_gibbs', myid
      return
      END

C  Added subroutine for global sync
C***********************************************************************
      Subroutine dntmc_sync(global)
C***********************************************************************
      Implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"

      logical global ! Input triggers a global sync event
      integer gidnow
C      integer world
C      integer g_sync, g_stop
C      integer ld(2), lds(2)
C      logical exit
C      integer status
C      integer synclist(2000)
C      integer i, test
C
      gidnow = ga_pgroup_get_default()

      call ga_sync()
C      if (global) then
C         world = ga_pgroup_get_world()
C         call ga_pgroup_set_default(world)
C      endif
C
C      if (.not.ga_create(MT_INT, ga_nnodes(), 1, "sync list",
C     &     1, -1, g_sync))
C     &   call errquit('dntmc_gibbs:g_sync create error', 0,
C     &                 GA_ERR)
C      if (.not.ga_create(MT_INT, ga_nnodes(), 1, "through list",
C     &     1, -1, g_stop))
C     &   call errquit('dntmc_gibbs:g_stop create error', 0,
C     &                 GA_ERR)
C      call ga_zero(g_sync)
C      call ga_zero(g_stop) 
C
C      ld(1) = 1
C      ld(2) = 1 
C      lds(1) = ga_nnodes()
C      lds(2) = 1
C      exit = .false.
C      do while (.not. exit)
C        status = 1
C        call ga_put(g_sync, ga_nodeid()+1, ga_nodeid()+1, 1, 1,
C     &              status, ld) 
C        call ga_get(g_sync, 1, ga_nnodes(), 1, 1, synclist, lds)
C        test = 0
C        do i = 1, ga_nnodes()
C          test = test + synclist(i)
C        end do
C        if (test .lt. ga_nnodes()) then
C          exit = .false.
C        else
C          exit = .true.
C          call ga_fill(g_stop, 1)
C        endif 
C      end do
C
C         test = 0
C      do i = 1, ga_nnodes()
C         test = test + synclist(i)
C      end do
C      if (test .eq. ga_nnodes()) then
C        exit = ga_destroy(g_sync)
C        exit = ga_destroy(g_stop)
C      endif
C
C      if (global) then
C        call ga_pgroup_set_default(gidnow)
C      endif
      write(6,*)'Made it past sync', ga_nodeid(), gidnow, global
      call flush(6)
      return
      end
C***********************************************************************
C  End of Subroutine dntmc_sync(global)
C***********************************************************************
C  End of added subroutine for sync
C  Added subroutines for restart
C***********************************************************************
      Subroutine dntmc_copy_filertdb(prcfil, rtdb, nmols, fprefix)
C***********************************************************************
      Implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"
#include "inp.fh"

      double precision uu, usum, uavr, usqrdsum, usqrd
      double precision rsim, rsimsum, rsimavr
      double precision seed, r, numbersum, numberavr
      integer number, iacceptt, iacceptv, iacceptr, iacceptpr, iacceptpi
      integer rsmcstps, i, j, ir
      integer natom, spc, rtdb
      double precision rsh
      integer nob, nmols
      character*4 prcfil
      character*256 fprefix
      character*256 fprefixcat
      integer ncent, k, rdim
      double precision coords 
      dimension coords(3, MXATOM*MXMOLS)

      dimension R(MXMOLS,MXATOM,3)
      dimension NUMBER(100)
      dimension NUMBERSUM(100)
      dimension NUMBERAVR(100)
      dimension SEED(1)

      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)

      nob = 100

      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(3(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.restart.',
     &      prcfil(1:4)

      OPEN(UNIT = 20,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'OLD')

C
C      OPEN(UNIT = 20, FILE = file_prefix(1:inp_strlen(file_prefix))
C     &    //'.restart.'//prcfil(1:4),
C     +    STATUS = "UNKNOWN")

         READ(20,*)SEED(1)
         READ(20,*)UU
         READ(20,*)USUM
         READ(20,*)UAVR
         READ(20,*)USQRDSUM
         READ(20,*)USQRD
         READ(20,*)RSIM
         READ(20,*)RSIMSUM
         READ(20,*)RSIMAVR
         READ(20,*)IACCEPTT
         READ(20,*)IACCEPTR
         READ(20,*)IACCEPTV
         READ(20,*)IACCEPTPR
         READ(20,*)IACCEPTPI
         READ(20,*)rsmcstps
         DO I = 1, NMOLS
            DO J = 1, NATOM(SPC(I))
               READ(20,*)(R(I,J,IR), IR = 1,3)
            ENDDO
         ENDDO
         
         DO I = 1, NOB
            READ(20,*)NUMBER(I)
            READ(20,*)NUMBERSUM(I)
            READ(20,*)RSH,NUMBERAVR(I)
         ENDDO

         CLOSE(20)
      
      if (.not. rtdb_put(rtdb, 'dntmc:seed:'//prcfil(1:4), MT_DBL, 1,
     &   SEED(1)))
     &   call errquit('dntmc_gibbs:  rtdb error for seed(1)',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:uu:'//prcfil(1:4), MT_DBL, 1, UU))
     &   call errquit('dntmc_gibbs:  rtdb error for UU',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:usum:'//prcfil(1:4), MT_DBL, 1,
     &   usum))
     &   call errquit('dntmc_gibbs:  rtdb error for usum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:uavr:'//prcfil(1:4), MT_DBL, 1,
     &   uavr))
     &   call errquit('dntmc_gibbs:  rtdb error for uavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:usqrdsum:'//prcfil(1:4), MT_DBL, 1
     &   ,usqrdsum))
     &   call errquit('dntmc_gibbs:  rtdb error for usqrdsum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:usqrd:'//prcfil(1:4), MT_DBL, 1,
     &   usqrd))
     &   call errquit('dntmc_gibbs:  rtdb error for usqrd',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsim:'//prcfil(1:4), MT_DBL, 1,
     &   rsim))
     &   call errquit('dntmc_gibbs:  rtdb error for rsim',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsimsum:'//prcfil(1:4), MT_DBL, 1,
     &   rsimsum))
     &   call errquit('dntmc_gibbs:  rtdb error for rsimsum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsimavr:'//prcfil(1:4), MT_DBL, 1,
     &   rsimavr))
     &   call errquit('dntmc_gibbs:  rtdb error for rsimavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptt:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptt))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptt',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptr:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptr))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptv:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptv))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptv',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptpr:'//prcfil(1:4), MT_INT,
     &   1, iacceptpr))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptpr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptpi:'//prcfil(1:4), MT_INT,
     &   1, iacceptpi))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptpi',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:number:'//prcfil(1:4), MT_INT, NOB
     &   ,number))
     &   call errquit('dntmc_gibbs:  rtdb error for number',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:numbersum:'//prcfil(1:4), MT_DBL,
     &   NOB,
     &   numbersum))
     &   call errquit('dntmc_gibbs:  rtdb error for numbersum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:numberavr:'//prcfil(1:4), MT_DBL,
     &   NOB,
     &   numberavr))
     &   call errquit('dntmc_gibbs:  rtdb error for numberavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsh:'//prcfil(1:4), MT_DBL, 1, rsh
     &   ))call errquit('dntmc_gibbs:  rtdb error for rsh',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsmcstps:'//prcfil(1:4), MT_INT, 1
     &   ,rsmcstps))
     &   call errquit('dntmc_gibbs:  rtdb error for rsmcstps',0,
     &   RTDB_ERR)

      ncent = 0
      Do i = 1, nmols
        Do j = 1, natom(spc(i))
          ncent = ncent + 1
          Do k = 1, 3
            coords(k, ncent) = r(i, j, k)
          end do
        end do
      end do

      rdim = 3*ncent

      if (.not. rtdb_put(rtdb, 'dntmc:r:'//prcfil(1:4), MT_DBL, rdim,
     &   coords))
     &   call errquit('dntmc_gibbs:  rtdb error for r',0,
     &   RTDB_ERR)

      return
      end

C***********************************************************************
C  End of Subroutine dntmc_copy_filertdb(prcfil, rtdb, nmols, fprefix)
C***********************************************************************

C***********************************************************************
      Subroutine dntmc_build_prcfil(prcfil, i)
C***********************************************************************
      Implicit none

      integer i
      character*4 prcfil

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

      return
      end

C***********************************************************************
C  End of Subroutine dntmc_build_prcfil(prcfil, i)
C***********************************************************************
 
C***********************************************************************
      Subroutine dntmc_write_prcfil(idmedium, groupnums, prcfil,ngroups)
C***********************************************************************
      Implicit none

#include "errquit.fh"

      integer i
      integer idmedium, groupnums, ngroups
      integer groupid
      character*4 prcfil

      dimension groupnums(ngroups)
     
      groupid = -1
      do i=1, ngroups
        if (groupnums(i) .eq. idmedium) then
          groupid = i
        endif
      enddo

      if (groupid .eq. -1)
     & call errquit('dntmc_gibbs:  failed to allocate idgroup',0,
     &   GA_ERR)
 
       call dntmc_build_prcfil(prcfil, groupid)

      return
      end

C***********************************************************************
C  End of Subroutine dntmc_write_prcfil(idmedium, groupnums, prcfil,ngroups)
C***********************************************************************

C***********************************************************************
      Subroutine dntmc_write_dntmc_nwchem(rtdb, nmols, prcfil)
C***********************************************************************
      Implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"

      integer ncent, i, j, k, nmols, natom, spc, rdim
      integer rtdb
      double precision r, coords, rcm, rs
      dimension coords(3, MXATOM*MXMOLS)
      character*4 prcfil

      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)

      ncent = 0
      Do i = 1, nmols
        Do j = 1, natom(spc(i))
          ncent = ncent + 1
          Do k = 1, 3
            coords(k, ncent) = r(i, j, k)
          end do
        end do
      end do

      rdim = 3*ncent

      if (.not. rtdb_put(rtdb, 'dntmc:r:'//prcfil(1:4), MT_DBL, rdim,
     &   coords))
     &   call errquit('dntmc_gibbs:  rtdb error for r',0,
     &   RTDB_ERR)

      return
      end
C***********************************************************************
C  End of Subroutine dntmc_write_dntmc_nwchem(rtdb, nmols, prcfil)
C***********************************************************************

C***********************************************************************
      Subroutine dntmc_read_dntmc_nwchem(rtdb, nmols, prcfil)
C***********************************************************************
      Implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"

      integer ncent, i, j, k, nmols, natom, spc, rdim
      integer rtdb
      double precision r, coords, rcm, rs
      dimension coords(3, MXATOM*MXMOLS)
      character*4 prcfil

      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)

      ncent = 0
      Do i = 1, nmols
        Do j = 1, natom(spc(i))
          ncent = ncent + 1
        end do
      end do

      rdim = 3*ncent

      if (.not. rtdb_get(rtdb, 'dntmc:r:'//prcfil(1:4), MT_DBL, rdim,
     &   coords))
     &   call errquit('dntmc_gibbs:  rtdb error for r',0,
     &   RTDB_ERR)

      ncent = 0
      Do i = 1, nmols
        Do j = 1, natom(spc(i))
          ncent = ncent + 1
          Do k = 1, 3
            r(i, j, k) = coords(k, ncent)
          end do
        end do
      end do

      return
      end
C***********************************************************************
C  End of Subroutine dntmc_read_dntmc_nwchem(rtdb, nmols, prcfil)
C***********************************************************************

C***********************************************************************
      Subroutine dntmc_readrestart(rtdb, uu, uavr, usqrd,
     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +  iacceptpi, rsh, prcfil, nmols, rsmcstps)
C***********************************************************************
      Implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"

      double precision uu, usum, uavr, usqrdsum, usqrd
      double precision rsim, rsimsum, rsimavr
      double precision seed, r, numbersum, numberavr
      double precision rmax, rmin, rcm, rs
      integer number, iacceptt, iacceptv, iacceptr, iacceptpr, iacceptpi
      integer rsmcstps
      double precision rsh
      integer rtdb
      integer nob, nmols
      character*4 prcfil

      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/AVR/USUM,USQRDSUM
      COMMON/NUMBER1/NUMBER(100)
      COMMON/NUMBER2/NUMBERSUM(100),NUMBERAVR(100)
      COMMON/SEEDS/SEED(1)

      nob = 100

      if (.not. rtdb_get(rtdb, 'dntmc:seed:'//prcfil(1:4), MT_DBL, 1,
     &   SEED(1)))
     &   call errquit('dntmc_gibbs:  rtdb error for seed(1)',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:uu:'//prcfil(1:4), MT_DBL, 1, UU))
     &   call errquit('dntmc_gibbs:  rtdb error for UU',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:usum:'//prcfil(1:4), MT_DBL, 1,
     &   usum))
     &   call errquit('dntmc_gibbs:  rtdb error for usum',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:uavr:'//prcfil(1:4), MT_DBL, 1,
     &   uavr))
     &   call errquit('dntmc_gibbs:  rtdb error for uavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:usqrdsum:'//prcfil(1:4), MT_DBL, 1
     &   ,usqrdsum))
     &   call errquit('dntmc_gibbs:  rtdb error for usqrdsum',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:usqrd:'//prcfil(1:4), MT_DBL, 1,
     &   usqrd))
     &   call errquit('dntmc_gibbs:  rtdb error for usqrd',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsim:'//prcfil(1:4), MT_DBL, 1,
     &   rsim))
     &   call errquit('dntmc_gibbs:  rtdb error for rsim',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsimsum:'//prcfil(1:4), MT_DBL, 1,
     &   rsimsum))
     &   call errquit('dntmc_gibbs:  rtdb error for rsimsum',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsimavr:'//prcfil(1:4), MT_DBL, 1,
     &   rsimavr))
     &   call errquit('dntmc_gibbs:  rtdb error for rsimavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:iacceptt:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptt))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptt',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:iacceptr:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptr))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptr',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:iacceptv:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptv))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptv',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:iacceptpr:'//prcfil(1:4), MT_INT,
     &   1, iacceptpr))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptpr',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:iacceptpi:'//prcfil(1:4), MT_INT,
     &   1, iacceptpi))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptpi',0,
     &   RTDB_ERR)

      call dntmc_read_dntmc_nwchem(rtdb, nmols, prcfil)

      if (.not. rtdb_get(rtdb, 'dntmc:number:'//prcfil(1:4), MT_INT, NOB
     &   ,number))
     &   call errquit('dntmc_gibbs:  rtdb error for number',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:numbersum:'//prcfil(1:4), MT_DBL,
     &   NOB,
     &   numbersum))
     &   call errquit('dntmc_gibbs:  rtdb error for numbersum',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:numberavr:'//prcfil(1:4), MT_DBL,
     &   NOB,
     &   numberavr))
     &   call errquit('dntmc_gibbs:  rtdb error for numberavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsh:'//prcfil(1:4), MT_DBL, 1, rsh
     &   ))call errquit('dntmc_gibbs:  rtdb error for rsh',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsmcstps:'//prcfil(1:4), MT_INT, 1
     &   ,rsmcstps))
     &   call errquit('dntmc_gibbs:  rtdb error for rsmcstps',0,
     &   RTDB_ERR)
 
      return
      end
C***********************************************************************
C  End of Subroutine dntmc_readrestart(rtdb, uu, uavr, usqrd,
C     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
C     +  iacceptpi, rsh, prcfil, nmols, rsmcstps)
C***********************************************************************

C***********************************************************************
      Subroutine dntmc_writerestart(rtdb, uu, uavr, usqrd,
     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +  iacceptpi, rsh, prcfil, nmols, rsmcstps)
C***********************************************************************
      Implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"

      double precision uu, usum, uavr, usqrdsum, usqrd
      double precision rsim, rsimsum, rsimavr
      double precision seed, r, numbersum, numberavr
      double precision rmin, rmax, rinc, rcm, rs
      integer number, iacceptt, iacceptv, iacceptr, iacceptpr, iacceptpi
      integer rsmcstps
      double precision rsh
      integer rtdb, i
      integer nob, nmols
      character*4 prcfil


      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/AVR/USUM,USQRDSUM
      COMMON/NUMBER1/NUMBER(100)
      COMMON/NUMBER2/NUMBERSUM(100),NUMBERAVR(100)
      COMMON/SEEDS/SEED(1)

      nob = 100

      if (.not. rtdb_put(rtdb, 'dntmc:seed:'//prcfil(1:4), MT_DBL, 1,
     &   SEED(1)))
     &   call errquit('dntmc_gibbs:  rtdb error for seed(1)',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:uu:'//prcfil(1:4), MT_DBL, 1, UU))
     &   call errquit('dntmc_gibbs:  rtdb error for UU',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:usum:'//prcfil(1:4), MT_DBL, 1,
     &   usum))
     &   call errquit('dntmc_gibbs:  rtdb error for usum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:uavr:'//prcfil(1:4), MT_DBL, 1,
     &   uavr))
     &   call errquit('dntmc_gibbs:  rtdb error for uavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:usqrdsum:'//prcfil(1:4), MT_DBL, 1
     &   ,usqrdsum))
     &   call errquit('dntmc_gibbs:  rtdb error for usqrdsum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:usqrd:'//prcfil(1:4), MT_DBL, 1,
     &   usqrd))
     &   call errquit('dntmc_gibbs:  rtdb error for usqrd',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsim:'//prcfil(1:4), MT_DBL, 1,
     &   rsim))
     &   call errquit('dntmc_gibbs:  rtdb error for rsim',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsimsum:'//prcfil(1:4), MT_DBL, 1,
     &   rsimsum))
     &   call errquit('dntmc_gibbs:  rtdb error for rsimsum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsimavr:'//prcfil(1:4), MT_DBL, 1,
     &   rsimavr))
     &   call errquit('dntmc_gibbs:  rtdb error for rsimavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptt:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptt))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptt',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptr:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptr))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptv:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptv))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptv',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptpr:'//prcfil(1:4), MT_INT,
     &   1, iacceptpr))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptpr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptpi:'//prcfil(1:4), MT_INT,
     &   1, iacceptpi))
     &   call errquit('dntmc_gibbs:  rtdb error for iacceptpi',0,
     &   RTDB_ERR)

      call dntmc_write_dntmc_nwchem(rtdb, nmols, prcfil)

      if (.not. rtdb_put(rtdb, 'dntmc:number:'//prcfil(1:4), MT_INT, NOB
     &   ,number))
     &   call errquit('dntmc_gibbs:  rtdb error for number',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:numbersum:'//prcfil(1:4), MT_DBL,
     &   NOB,
     &   numbersum))
     &   call errquit('dntmc_gibbs:  rtdb error for numbersum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:numberavr:'//prcfil(1:4), MT_DBL,
     &   NOB,
     &   numberavr))
     &   call errquit('dntmc_gibbs:  rtdb error for numberavr',0,
     &   RTDB_ERR)

         RINC = (RMAX - RMIN)/DBLE(NOB)
         RSH = RMIN

         DO I = 1, NOB
            RSH = RSH + RINC
         ENDDO

      if (.not. rtdb_put(rtdb, 'dntmc:rsh:'//prcfil(1:4), MT_DBL, 1, rsh
     &   ))call errquit('dntmc_gibbs:  rtdb error for rsh',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsmcstps:'//prcfil(1:4), MT_INT, 1
     &   ,rsmcstps))
     &   call errquit('dntmc_gibbs:  rtdb error for rsmcstps',0,
     &   RTDB_ERR)

      return
      end

C***********************************************************************
C  End of Subroutine dntmc_writerestart(rtdb, uu, uavr, usqrd,
C     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
C     +  iacceptpi, rsh, prcfil, nmols, rsmcstps)
C***********************************************************************


C************************************************************************
      Subroutine dntmc_writefiles(uu, uavr, usqrd,
     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, fprefix)
C***********************************************************************
      Implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "inp.fh"

      double precision uu, usum, uavr, usqrdsum, usqrd
      double precision rsim, rsimsum, rsimavr
      double precision seed, r, numbersum, numberavr
      double precision rmin, rmax, rinc, rcm, rs
      integer number, iacceptt, iacceptv, iacceptr, iacceptpr, iacceptpi
      integer rsmcstps, i, j, ir
      integer natom, spc
      double precision rsh
      integer nob, nmols
      character*4 prcfil
      character*256 fprefix
      character*256 fprefixcat

      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/AVR/USUM,USQRDSUM
      COMMON/NUMBER1/NUMBER(100)
      COMMON/NUMBER2/NUMBERSUM(100),NUMBERAVR(100)
      COMMON/SEEDS/SEED(1)
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)

      nob = 100

      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(3(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.restart.',
     &      prcfil(1:4)

      OPEN(UNIT = 20,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')

C
C      OPEN(UNIT = 20, FILE = file_prefix(1:inp_strlen(file_prefix))
C     &    //'.restart.'//prcfil(1:4),
C     +    STATUS = "UNKNOWN")

         WRITE(20,*)SEED(1)
         WRITE(20,*)UU
         WRITE(20,*)USUM
         WRITE(20,*)UAVR
         WRITE(20,*)USQRDSUM
         WRITE(20,*)USQRD
         WRITE(20,*)RSIM
         WRITE(20,*)RSIMSUM
         WRITE(20,*)RSIMAVR
         WRITE(20,*)IACCEPTT
         WRITE(20,*)IACCEPTR
         WRITE(20,*)IACCEPTV
         WRITE(20,*)IACCEPTPR
         WRITE(20,*)IACCEPTPI
         WRITE(20,*)rsmcstps
         DO I = 1, NMOLS
            DO J = 1, NATOM(SPC(I))
               WRITE(20,*)(R(I,J,IR), IR = 1,3)
            ENDDO
         ENDDO

         NOB = 100

         RINC = (RMAX - RMIN)/DBLE(NOB)
         RSH = RMIN

         DO I = 1, NOB
            RSH = RSH + RINC
            WRITE(20,*)NUMBER(I)
            WRITE(20,*)NUMBERSUM(I)
            WRITE(20,*)RSH,NUMBERAVR(I)
         ENDDO

         CLOSE(20)

      return
      end

C***********************************************************************
C  End of Subroutine dntmc_writefiles(uu, uavr, usqrd,
C     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
C     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, fprefix)
C*********************************************************************** 

C  End of Subroutines for restart

C***********************************************************************
      SUBROUTINE dntmc_MOVE_RWEIGHT(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,
     +               ACCEPTPR,ACCEPTPI,UOLD, rtdb, nmols)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"

      double precision temp  ! input
      LOGICAL ACCEPTT        ! output
      LOGICAL ACCEPTR        ! output
      LOGICAL ACCEPTV        ! output
      LOGICAL ACCEPTPR       ! output
      LOGICAL ACCEPTPI       ! output
      double precision uold  ! input/output
      integer rtdb           ! input
      integer nmols          ! input

      logical task_energy
      external task_energy

      double precision r, rcm, rs
      double precision rmin, rsim, rmax
      double precision rclcm
      double precision e, k, esqrd, kb, pi
      INTEGER NATOM
      double precision disp, alpha
      INTEGER ROTF
      INTEGER SPC
      INTEGER ANDX
      INTEGER NMIN,NCURR,NMAX
      integer l
      double precision seed
      double precision rnnum
      integer idum

      double precision delta
      double precision comp2, comp3, comp2r, comp3r, press
      double precision rsimsq, radiuss, deltar, rnew
      double precision voli, volf, deltavol, rnewsq
      double precision rx, ry, rz, dist, arg, tmp, rann
      double precision unew, deltau, beta, sdbeta, cdbeta
      integer i, j, ir, ia1, ia2, ia3

      double precision dntmc_ran3
C Added array for configuration files
      double precision DISTcalc(MXMOLS)
C End added array
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/COORD2/RCLCM(3)
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/TRPS/DISP,ALPHA,ROTF
      COMMON/SPECIES/SPC(MXMOLS)          
      COMMON/ANDX/ANDX(MXSP,MXSP)      
      COMMON/GRAND/NMIN,NCURR,NMAX
      COMMON/MSTP/L
      COMMON/SEEDS/SEED(1)
      COMMON/RNNUMS/RNNUM(1)
      COMMON/SEED/IDUM

      DIMENSION DELTA(3)

      COMP2 = 0.D+00
      COMP3 = 0.D+00
      COMP2R = 0.D+00
      COMP3R = 0.D+00

c Do these need to more dependant on the actual temperature 
c of the simulation?
      PRESS = 1.7D-05   !kcal/mole/A^3  P = 888 torr  taken from UMR data

      RSIMSQ = RSIM*RSIM

      ACCEPTV = .FALSE.

C Save original coordinates

      DO I = 1, NMOLS
           DO J = 1, NATOM(SPC(I))
              RS(I,J,1) = R(I,J,1)
              RS(I,J,2) = R(I,J,2)
              RS(I,J,3) = R(I,J,3)
           END DO
      END DO

C Lets do translation

      DO I = 1, NMOLS
C     PICK DISPLACEMENTS ALONG X,Y,Z
         DELTA(1) = DISP*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
         DELTA(2) = DISP*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
         DELTA(3) = DISP*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) + DELTA(1)
            R(I,J,2) = R(I,J,2) + DELTA(2)
            R(I,J,3) = R(I,J,3) + DELTA(3)
         END DO
      END DO

C Center Molecule

      CALL dntmc_CLUSTERCMASS(NMOLS)

      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

C  Get new rsim

      DO I = 1, NMOLS
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DISTcalc(I) = RX*RX + RY*RY + RZ*RZ
      END DO
      DIST = DISTcalc(1)
      DO I = 2, NMOLS
         IF (DISTcalc(I) .gt. DIST) then
            DIST = DISTcalc(I)
         endif
      enddo

      rnew = DSQRT(DIST)

C Reject based on volume

      if ((rnew .gt. rmax).or.(rnew .lt. rmin)) goto 20

C Do Metropolis on translatonal move

      call dntmc_setcoords(rtdb, r, nmols, natom, spc)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          write(12,*) 'Energy did not converge'
          unew = -0.5                        ! set the energy arbitrarily high - H atom
        else
          write(12,*) 'Energy converged on the second try'
          if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $      call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
        endif
      else
        if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $     call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
      endif
c      if (.not. task_energy(rtdb)) call errquit
c     $     ('dntmc_move: energy failed',0, GEOM_ERR)
c      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
c     $     call errquit('dntmc_move: rtdb?',0, RTDB_ERR)

      write(12,*) 'coordinates with an energy evaluation'
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            WRITE(12,*)(R(I,J,IR),IR=1,3)
C Added write for configuration files
            WRITE(10,*)(R(I,J,IR),IR=1,3)
C End Added write
         END DO
      END DO
c
c convert to kcal/mol
c
C Added write for configuration files
      WRITE(10,*)'Ucalc = ', Unew,' hartree'
      WRITE(10,*)'Rsim = ', RNEW,' Angstrom'

C End added write 
      unew = unew*627.5095D+00
      write(12,*) 'unew = ',unew

c      WRITE(12,*)'after translation from psum'
c      WRITE(12,*)'UNEW =',UNEW
c      WRITE(12,*)'UOLD =',UOLD

      ARG = 0.D+00

      DELTAU = UNEW - UOLD

      ARG = DELTAU/(kb*TEMP)

      TMP = ((rnew/rsim)**2.D+00)*DEXP(-ARG)

c      CALL dntmc_RAN(1)

      RANN = dntmc_RAN3(IDUM)

c      WRITE(12,*)'TMP =',TMP
c      WRITE(12,*)'RANN =',RANN
      
      IF (TMP.GE.1.D+00) GOTO 15

      IF (TMP.LT.RANN) GOTO 20


C     ACCEPT

  15  CONTINUE

      ACCEPTT = .TRUE.
      UOLD = UNEW
      rsim = rnew

C  Save Coordinates

      DO I = 1, NMOLS
           DO J = 1, NATOM(SPC(I))
              RS(I,J,1) = R(I,J,1)
              RS(I,J,2) = R(I,J,2)
              RS(I,J,3) = R(I,J,3)
           END DO
      END DO

      GOTO 25

C  Reject

  20  CONTINUE

      ACCEPTT = .FALSE.

C Reload Coordinates

      DO I = 1, NMOLS
           DO J = 1, NATOM(SPC(I))
              R(I,J,1) = RS(I,J,1)
              R(I,J,2) = RS(I,J,2)
              R(I,J,3) = RS(I,J,3)
           END DO
      END DO

   25  CONTINUE
C  Lets do rotation

      DO I = 1, NMOLS

C     PICK AN ANGLE
      BETA = ALPHA*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)

      SDBETA = DSIN(BETA)
      CDBETA = DCOS(BETA)

C     PICK AN AXIS ABOUT WHICH TO ROTATE MOLECULE I
      IA1 = INT(3.D+00*dntmc_RAN3(IDUM) + 1)
      IA2 = IA1 + 1
      IA3 = IA1 + 2
      IF(IA2.GT.3) IA2 = IA2 - 3
      IF(IA3.GT.3) IA3 = IA3 - 3

C     ROTATE THE MOLECULE

      CALL dntmc_CMASS(I)

         DO J = 1, NATOM(SPC(I))
            COMP2 = (R(I,J,IA2) - RCM(I,IA2))
            COMP3 = (R(I,J,IA3) - RCM(I,IA3))

            COMP2R=CDBETA*COMP2+SDBETA*COMP3
            COMP3R=-SDBETA*COMP2+CDBETA*COMP3

            R(I,J,IA2)= RCM(I,IA2) + COMP2R
            R(I,J,IA3)= RCM(I,IA3) + COMP3R
         END DO
      END DO

C Center Molecule

      CALL dntmc_CLUSTERCMASS(NMOLS)

      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

C  Get new rsim

      DO I = 1, NMOLS
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DISTcalc(I) = RX*RX + RY*RY + RZ*RZ
      END DO
      DIST = DISTcalc(1)
      DO I = 2, NMOLS
         IF (DISTcalc(I) .gt. DIST) then
            DIST = DISTcalc(I)
         endif
      enddo

      rnew = DSQRT(DIST)

C Reject based on volume

      if ((rnew .gt. rmax).or.(rnew .lt. rmin)) goto 35

C Do Metropolis on rotation

      call dntmc_setcoords(rtdb, r, nmols, natom, spc)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          write(12,*) 'Energy did not converge'
          unew = -0.5                        ! set the energy arbitrarily high - H atom
        else
          write(12,*) 'Energy converged on the second try'
          if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $      call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
        endif
      else
        if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $     call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
      endif
c      if (.not. task_energy(rtdb)) call errquit
c     $     ('dntmc_move: energy failed',0, GEOM_ERR)
c      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
c     $     call errquit('dntmc_move: rtdb?',0, RTDB_ERR)

      write(12,*) 'coordinates with an energy evaluation'
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            WRITE(12,*)(R(I,J,IR),IR=1,3)
C Added write for configuration files
            WRITE(10,*)(R(I,J,IR),IR=1,3)
C End Added write
         END DO
      END DO
c
c convert to kcal/mol
c
C Added write for configuration files
      WRITE(10,*)'Ucalc = ', Unew,' hartree'
      WRITE(10,*)'Rsim = ', RNEW,' Angstrom'
C End added write 
      unew = unew*627.5095D+00
      write(12,*) 'unew = ',unew

c      WRITE(12,*)'after rotation from psum'
c      WRITE(12,*)'UNEW =',UNEW
c      WRITE(12,*)'UOLD =',UOLD

      DELTAU = UNEW - UOLD

      ARG = DELTAU/(kb*TEMP)

      TMP = ((rnew/rsim)**2.D+00)*DEXP(-ARG)


      RANN = dntmc_RAN3(IDUM)

      IF (TMP.GE.1.D+00) GOTO 30 

      IF (TMP.LT.RANN) GOTO 35

C Accept
   30 CONTINUE

      ACCEPTR = .TRUE.
      UOLD = UNEW
      RSIM = RNEW

C  Save Coordinates

      DO I = 1, NMOLS
           DO J = 1, NATOM(SPC(I))
              RS(I,J,1) = R(I,J,1)
              RS(I,J,2) = R(I,J,2)
              RS(I,J,3) = R(I,J,3)
           END DO
      END DO

      GOTO 40

C  Reject
   35 CONTINUE

      ACCEPTR = .FALSE.

C  Reload Coordinates

      DO I = 1, NMOLS
           DO J = 1, NATOM(SPC(I))
              R(I,J,1) = RS(I,J,1)
              R(I,J,2) = RS(I,J,2)
              R(I,J,3) = RS(I,J,3)
           END DO
      END DO

   40 CONTINUE
C  Finish up

C Added write for configuration files
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            WRITE(11,*)(R(I,J,IR),IR=1,3)
         enddo
      enddo
         WRITE(11,*)'Ucalc = ',(Uold/627.5095D+00),' hartree'
         WRITE(11,*)'Rsim = ',rsim,' Angstrom'
C  End added write
      RETURN
      END
C***********************************************************************
C  End of SUBROUTINE dntmc_MOVE_RWEIGHT(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,
C     +               ACCEPTPR,ACCEPTPI,UOLD, rtdb, nmols)
C***********************************************************************


C***********************************************************************
      SUBROUTINE dntmc_MOVE(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
     +               ,ACCEPTPI,UOLD, rtdb, nmols)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"

      double precision temp  ! input
      LOGICAL ACCEPTT        ! output
      LOGICAL ACCEPTR        ! output
      LOGICAL ACCEPTV        ! output
      LOGICAL ACCEPTPR       ! output
      LOGICAL ACCEPTPI       ! output
      double precision uold  ! input/output
      integer rtdb           ! input
      integer nmols          ! input

      logical task_energy
      external task_energy

      double precision r, rcm, rs
      double precision rmin, rsim, rmax
      double precision rclcm
      double precision e, k, esqrd, kb, pi
      INTEGER NATOM
      double precision disp, alpha
      INTEGER ROTF
      INTEGER SPC
      INTEGER ANDX
      INTEGER NMIN,NCURR,NMAX
      integer l
      double precision seed
      double precision rnnum
      integer idum

      double precision delta
      double precision comp2, comp3, comp2r, comp3r, rstep, press
      double precision rsimsq, radiuss, deltar, rnew
      double precision voli, volf, deltavol, rnewsq
      double precision rx, ry, rz, dist, arg, tmp, rann
      double precision unew, deltau, beta, sdbeta, cdbeta
      integer i, j, ir, ia1, ia2, ia3

      double precision dntmc_ran3
C Added array for configuration files
      double precision DISTcalc(MXMOLS)
C End added array
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/COORD2/RCLCM(3)
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/TRPS/DISP,ALPHA,ROTF
      COMMON/SPECIES/SPC(MXMOLS)          
      COMMON/ANDX/ANDX(MXSP,MXSP)      
      COMMON/GRAND/NMIN,NCURR,NMAX
      COMMON/MSTP/L
      COMMON/SEEDS/SEED(1)
      COMMON/RNNUMS/RNNUM(1)
      COMMON/SEED/IDUM

      DIMENSION DELTA(3)

      COMP2 = 0.D+00
      COMP3 = 0.D+00
      COMP2R = 0.D+00
      COMP3R = 0.D+00
      RSTEP = (RMAX - RMIN)

c Do these need to more dependant on the actual temperature 
c of the simulation?
      PRESS = 1.7D-05   !kcal/mole/A^3  P = 888 torr  taken from UMR data

      RSIMSQ = RSIM*RSIM

C*****************************************************
C     CHANGE THE SIMULATION VOLUME
C*****************************************************

      RADIUSS = RSIM

  9   CONTINUE

c      CALL dntmc_RAN(1)
C change in the radius
      DELTAR = RSTEP*(dntmc_RAN3(IDUM) - 0.5D+00)

c     WRITE(12,*)'DELTAR =',DELTAR,'at step L =',L 

      RNEW = RSIM + DELTAR

      IF(RNEW.GT.RMAX)THEN
        RNEW = 2.D+00*RMAX - RNEW
      ELSEIF(RNEW.LT.RMIN)THEN
        RNEW = 2.D+00*RMIN - RNEW
      ENDIF

      VOLI = (RADIUSS**3.D+00)*(4.D+00*PI)/(3.D+00)

      VOLF = (RNEW**3.D+00)*(4.D+00*PI)/(3.D+00)

      DELTAVOL = VOLF - VOLI

      RNEWSQ = RNEW*RNEW

C find center of mass for each molecule in the current config
      DO I = 1, NCURR
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DIST = RX*RX + RY*RY + RZ*RZ

C reject based on volume
         IF(DIST.GT.RNEWSQ) goto 11

      END DO

C******************************************************
C******************************************************
C     Do metropolis on volume change only
C******************************************************
C******************************************************
      
      ARG = PRESS*(DELTAVOL)/(kb*TEMP)

c     WRITE(12,*)'PRESSURE =',PRESS
c     WRITE(12,*)'DELTAVOL =',DELTAVOL
c     WRITE(12,*)'ARG =',ARG

      TMP = ((VOLF/VOLI)**(2.D+00/3.D+00))*DEXP(-ARG)

c      CALL dntmc_RAN(1)

      RANN = dntmc_RAN3(IDUM)

c     WRITE(12,*)'TMP =',TMP
c     WRITE(12,*)'RANN =',RANN

      IF (TMP.GE.1.D+00) GOTO 10 

      IF (TMP.LT.RANN) GOTO 11

C**********************************************
C     ACCEPT
C**********************************************

  10  CONTINUE

      ACCEPTV = .TRUE.
      
c      WRITE(12,*)'VOLUME ACCEPTED'

      RSIM = RNEW

      VOLI = VOLF

c     WRITE(12,*)'RSIM =',RSIM,'at step L =',L   

      GOTO 12

C**********************************************
C     REJECT
C**********************************************

  11  CONTINUE

      ACCEPTV = .FALSE.

c      WRITE(12,*)'VOLUME REJECTED'

      RSIM = RADIUSS

      VOLF = VOLI
    

  12  CONTINUE

      RSIMSQ = RSIM*RSIM

C*****************************************************
C     SAVE THE COORDINATES OF ALL THE MOLECULES
C*****************************************************
  
      DO I = 1, NCURR 
           DO J = 1, NATOM(SPC(I))
              RS(I,J,1) = R(I,J,1)
              RS(I,J,2) = R(I,J,2)
              RS(I,J,3) = R(I,J,3)
           END DO
      END DO      

C******************************************************
C     TRANSLATE ALL MOLECULES WITHIN SIMULATION VOLUME
C******************************************************

      GOTO 14

C should not ever get here
  13  CONTINUE

C************************************************************
C     RELOAD COORDINATES OF ALL THE MOLECULES 
C************************************************************

c      WRITE(12,*)'in move routine'
c      WRITE(12,*)'coordinates are:'

      DO I = 1, NCURR
         DO J = 1, NATOM(SPC(I))
             R(I,J,1) = RS(I,J,1)
             R(I,J,2) = RS(I,J,2)
             R(I,J,3) = RS(I,J,3)
c            WRITE(12,*)(R(I,J,IR),IR = 1,3)
         END DO   
      END DO

  14  CONTINUE

      DO I = 1, NCURR
C     PICK DISPLACEMENTS ALONG X,Y,Z
c      CALL dntmc_RAN(1)
         DELTA(1) = DISP*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
c      CALL dntmc_RAN(1)
         DELTA(2) = DISP*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
c      CALL dntmc_RAN(1)
         DELTA(3) = DISP*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) + DELTA(1)
            R(I,J,2) = R(I,J,2) + DELTA(2)
            R(I,J,3) = R(I,J,3) + DELTA(3)
         END DO
      END DO

      CALL dntmc_CLUSTERCMASS(NCURR)

      DO I = 1, NCURR
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

      DO I = 1, NCURR
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DIST = RX*RX + RY*RY + RZ*RZ

         IF(DIST.GT.RSIMSQ) goto 16

      END DO

C******************************************************
C******************************************************
C     Do metropolis on translational part only
C******************************************************
C******************************************************

c      WRITE(12,*)'translation before psum' 

c     CALL dntmc_PSUM(NCURR,UNEW) 
      call dntmc_setcoords(rtdb, r, nmols, natom, spc)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          write(12,*) 'Energy did not converge'
          unew = -0.5                        ! set the energy arbitrarily high - H atom
        else
          write(12,*) 'Energy converged on the second try'
          if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $      call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
        endif
      else 
        if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $     call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
      endif
c      if (.not. task_energy(rtdb)) call errquit
c     $     ('dntmc_move: energy failed',0, GEOM_ERR)
c      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
c     $     call errquit('dntmc_move: rtdb?',0, RTDB_ERR)

      write(12,*) 'coordinates with an energy evaluation'
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            WRITE(12,*)(R(I,J,IR),IR=1,3)
C Added write for configuration files
            WRITE(10,*)(R(I,J,IR),IR=1,3)
C End Added write
         END DO
      END DO
c
c convert to kcal/mol
c
C Added write for configuration files
      WRITE(10,*)'Ucalc = ', Unew,' hartree'
      WRITE(10,*)'Rsim = ', RSIM,' Angstrom'
      DO I = 1, NMOLS
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DISTcalc(I) = RX*RX + RY*RY + RZ*RZ
      END DO 
      DIST = DISTcalc(1)
      DO I = 2, NMOLS
         IF (DISTcalc(I) .gt. DIST) then
            DIST = DISTcalc(I)
         endif
      enddo
      WRITE(10,*)'Rsim (Calculated) = ',DSQRT(DIST),' Angstrom'
C End added write 
      unew = unew*627.5095D+00
      write(12,*) 'unew = ',unew

c      WRITE(12,*)'after translation from psum'
c      WRITE(12,*)'UNEW =',UNEW
c      WRITE(12,*)'UOLD =',UOLD

      ARG = 0.D+00

      DELTAU = UNEW - UOLD

      ARG = DELTAU/(kb*TEMP)
      
      TMP = DEXP(-ARG)

c      CALL dntmc_RAN(1)

      RANN = dntmc_RAN3(IDUM)

c      WRITE(12,*)'TMP =',TMP
c      WRITE(12,*)'RANN =',RANN

      IF (TMP.GE.1.D+00) GOTO 15 

      IF (TMP.LT.RANN) GOTO 16


C**********************************************
C     ACCEPT
C**********************************************

  15  CONTINUE
  
      ACCEPTT = .TRUE.
      
c      WRITE(12,*)'TRANSLATION ACCEPTED'
  
      UOLD = UNEW
    
      GOTO 17

C**********************************************
C     REJECT
C**********************************************
  16  CONTINUE

      ACCEPTT = .FALSE.

c      WRITE(12,*)'TRANSLATION REJECTED coords are:'
    
C     RELOAD OLD COORDINATES AFTER TRANSLATION
      
      DO I = 1, NCURR 
           DO J = 1, NATOM(SPC(I))
              R(I,J,1) = RS(I,J,1)
              R(I,J,2) = RS(I,J,2)
              R(I,J,3) = RS(I,J,3)
c               WRITE(12,*)(R(I,J,IR),IR = 1,3)
           END DO
      END DO           

  17  CONTINUE

C*****************************************************
C     SAVE THE COORDINATES OF ALL THE MOLECULES
C*****************************************************
  
      DO I = 1, NCURR 
           DO J = 1, NATOM(SPC(I))
              RS(I,J,1) = R(I,J,1)
              RS(I,J,2) = R(I,J,2)
              RS(I,J,3) = R(I,J,3)
           END DO
      END DO

C******************************************************
C     ROTATE ALL MOLECULES WITHIN SIMULATION VOLUME
C******************************************************

      GOTO 19

  18  CONTINUE

C************************************************************  
C     RELOAD COORDINATES OF ALL THE MOLECULES
C************************************************************

      DO I = 1, NCURR  
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = RS(I,J,1)
            R(I,J,2) = RS(I,J,2)
            R(I,J,3) = RS(I,J,3)
         END DO   
      END DO

  19  CONTINUE

      DO I = 1, NCURR

C     PICK AN ANGLE
c      CALL dntmc_RAN(1)
      BETA = ALPHA*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)

c      WRITE(12,*)'ALPHA =',ALPHA
c      WRITE(12,*)'ROTATION ANGLE =',BETA,' radians'

      SDBETA = DSIN(BETA)
      CDBETA = DCOS(BETA)

C     PICK AN AXIS ABOUT WHICH TO ROTATE MOLECULE I
c      CALL dntmc_RAN(1)
      IA1 = INT(3.D+00*dntmc_RAN3(IDUM) + 1)   
      IA2 = IA1 + 1
      IA3 = IA1 + 2
      IF(IA2.GT.3) IA2 = IA2 - 3
      IF(IA3.GT.3) IA3 = IA3 - 3

c     PRINT*,'ROT.AXIS(xyz->123) =',IA1,IA2,IA3
C     ROTATE THE MOLECULE
      
      CALL dntmc_CMASS(I)

         DO J = 1, NATOM(SPC(I))
            COMP2 = (R(I,J,IA2) - RCM(I,IA2))
            COMP3 = (R(I,J,IA3) - RCM(I,IA3))

            COMP2R=CDBETA*COMP2+SDBETA*COMP3
            COMP3R=-SDBETA*COMP2+CDBETA*COMP3

            R(I,J,IA2)= RCM(I,IA2) + COMP2R
            R(I,J,IA3)= RCM(I,IA3) + COMP3R
         END DO
      
      END DO

      CALL dntmc_CLUSTERCMASS(NCURR)

      DO I = 1, NCURR
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

      DO I = 1, NCURR
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DIST = RX*RX + RY*RY + RZ*RZ

         IF(DIST.GT.RSIMSQ) goto 18
      END DO

C******************************************************
C******************************************************
C     Do metropolis on rotational part only
C******************************************************
C******************************************************
      
c      CALL dntmc_PSUM(NCURR,UNEW)     
      call dntmc_setcoords(rtdb, r, nmols, natom, spc)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          write(12,*) 'Energy did not converge'
          unew = -0.5                        ! set the energy arbitrarily high - H atom
        else
          write(12,*) 'Energy converged on the second try'
          if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $      call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
        endif
      else 
        if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $     call errquit('dntmc_gibbs: rtdb?',0, RTDB_ERR)
      endif
c      if (.not. task_energy(rtdb)) call errquit
c     $     ('dntmc_move: energy failed',0, GEOM_ERR)
c      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
c     $     call errquit('dntmc_move: rtdb?',0, RTDB_ERR)

      write(12,*) 'coordinates with an energy evaluation'
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            WRITE(12,*)(R(I,J,IR),IR=1,3)
C Added write for configuration files
            WRITE(10,*)(R(I,J,IR),IR=1,3)
C End Added write
         END DO
      END DO
c
c convert to kcal/mol
c
C Added write for configuration files
      WRITE(10,*)'Ucalc = ', Unew,' hartree'
      WRITE(10,*)'Rsim = ', RSIM,' Angstrom'
      DO I = 1, NMOLS
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DISTcalc(I) = RX*RX + RY*RY + RZ*RZ
      END DO 
      DIST = DISTcalc(1)
      DO I = 2, NMOLS
         IF (DISTcalc(I) .gt. DIST) then
            DIST = DISTcalc(I)
         endif
      enddo
      WRITE(10,*)'Rsim (Calculated) = ',DSQRT(DIST),' Angstrom'
C End added write 
      unew = unew*627.5095D+00
      write(12,*) 'unew = ',unew

c      WRITE(12,*)'after rotation from psum'
c      WRITE(12,*)'UNEW =',UNEW
c      WRITE(12,*)'UOLD =',UOLD

      DELTAU = UNEW - UOLD

      ARG = DELTAU/(kb*TEMP)
      
      TMP = DEXP(-ARG)

c      CALL dntmc_RAN(1)

      RANN = dntmc_RAN3(IDUM)

c      WRITE(12,*)'TMP =',TMP
c      WRITE(12,*)'RANN =',RANN

      IF (TMP.GE.1.D+00) GOTO 20

      IF (TMP.LT.RANN) GOTO 21

C**********************************************
C     ACCEPT
C**********************************************

  20  CONTINUE

      ACCEPTR = .TRUE.
      
c      WRITE(12,*)'ROTATION ACCEPTED'

      UOLD = UNEW
    
      GOTO 22

C**********************************************
C     REJECT
C**********************************************

  21  CONTINUE   

      ACCEPTR = .FALSE.

c      WRITE(12,*)'ROTATION REJECTED coords are:'

C     RELOAD OLD COORDINATES
      
      DO I = 1, NCURR 
           DO J = 1, NATOM(SPC(I))
              R(I,J,1) = RS(I,J,1)
              R(I,J,2) = RS(I,J,2)
              R(I,J,3) = RS(I,J,3)
c               WRITE(12,*)(R(I,J,IR),IR=1,3)
           END DO
      END DO         

  22  CONTINUE

C Added write for configuration files
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            WRITE(11,*)(R(I,J,IR),IR=1,3)
         enddo
      enddo
      DO I = 1, NMOLS
         CALL dntmc_CMASS(I)
         RX = RCM(I,1)
         RY = RCM(I,2)
         RZ = RCM(I,3)
         DISTcalc(I) = RX*RX + RY*RY + RZ*RZ
      END DO
      DIST = DISTcalc(1)
      DO I = 2, NMOLS
         IF (DISTcalc(I) .gt. DIST) then
            DIST = DISTcalc(I)
         endif
      enddo
         WRITE(11,*)'Ucalc = ',(Uold/627.5095D+00),' hartree'
         WRITE(11,*)'Rsim = ',rsim,' Angstrom'
         WRITE(11,*)'Rsim (Calculated) = ',DSQRT(DIST),' Angstrom'
C  End added write
      RETURN
      END
      
C***********************************************************************
      SUBROUTINE dntmc_PSUM(NMOLS,UU)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "errquit.fh"
#include "gibbsP.fh"

      integer nmols        ! input
      double precision uu  ! output

      double precision r, rcm, rs
      integer natom
      integer spc
      integer andx
      integer lndx,  lndxx
      double precision LAMBDA
      integer ab
      integer l
      double precision e, k, esqrd, kb, pi
      integer jindex
      double precision r3

      double precision LAM
      double precision MU
c      real*8 MNUX,MNUY,MNUZ,MTX,MTY,MTZ,MTMAG
      double precision NUX,NUY,NUZ,MUX,MUY,MUZ,MX,MY,MZ,MUTOL
      INTEGER NSCF,SCF,PLRZ,PLRZA
      double precision v, ef, efwx, efwy, efwz
      double precision efix, efiy, efiz
      double precision efx, efy, efz
      double precision pol, count, vh2o, elam
      double precision upol1, rxiijj, ryiijj, rziijj, riijjsq, riijj
      double precision efux, efuy, efuz, efldmag
      double precision sum, dotx, doty, dotz, dot, dmu, x, udot1
      integer i, j, ii, jj

      double precision dntmc_efmag
      double precision dntmc_VLjones
      double precision dntmc_Vcoul

      PARAMETER (NSCF   = 200)            

      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)
      COMMON/ANDX/ANDX(MXSP,MXSP)
      COMMON/LAMBNDX/LNDX(MXMOLS),LAMBDA,LNDXX(9,9)
      COMMON/AORB/AB
      COMMON/MSTP/L      
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/FORMAT/JINDEX
      COMMON/CUBE/R3(MXMOLS,MXMOLS)      

      DIMENSION V(MXMOLS,MXMOLS,MXATOM,MXATOM)
      DIMENSION EF(MXMOLS,MXMOLS,MXATOM,MXATOM,3)
      DIMENSION EFWX(MXMOLS),EFWY(MXMOLS),EFWZ(MXMOLS)
      DIMENSION EFIX(MXMOLS),EFIY(MXMOLS),EFIZ(MXMOLS)
      DIMENSION EFX(MXMOLS),EFY(MXMOLS),EFZ(MXMOLS)
      DIMENSION NUX(MXMOLS,MXMOLS),NUY(MXMOLS,MXMOLS)
      DIMENSION NUZ(MXMOLS,MXMOLS)
      DIMENSION MUX(NSCF,MXMOLS),MUY(NSCF,MXMOLS)
      DIMENSION MUZ(NSCF,MXMOLS),MU(NSCF,MXMOLS)
      DIMENSION MX(MXMOLS),MY(MXMOLS),MZ(MXMOLS)
      DIMENSION POL(MXSP),COUNT(MXMOLS)            
      DIMENSION Vh2o(MXATOM,MXATOM)
      DIMENSION LAM(2,0:1,0:1,0:1)
      DIMENSION ELAM(0:1,0:1)


      LAM(1,0,0,0) = 1.D+00
      LAM(1,0,0,1) = 1.D+00
      LAM(1,0,1,0) = 1.D+00
      LAM(1,0,1,1) = 1.D+00 
      LAM(1,1,0,0) = 1.D+00
      LAM(1,1,0,1) = 1.D+00
      LAM(1,1,1,0) = 1.D+00
      LAM(1,1,1,1) = 1.D+00
      LAM(2,0,0,0) = 1.D+00
      LAM(2,0,0,1) = 1.D+00
      LAM(2,0,1,0) = LAMBDA
      LAM(2,0,1,1) = LAMBDA
      LAM(2,1,0,0) = LAMBDA
      LAM(2,1,0,1) = LAMBDA
      LAM(2,1,1,0) = LAMBDA
      LAM(2,1,1,1) = 1.D+00

      ELAM(0,0) = 1.D+00
      ELAM(0,1) = 1.D+00
      ELAM(1,0) = 1.D+00
      ELAM(1,1) = 1.D+00     
      
      UU = 0.D+00
      Upol1 = 0.D+00

      DO I = 1, NMOLS
         EFWX(I) = 0.D+00
         EFWY(I) = 0.D+00
         EFWZ(I) = 0.D+00
      ENDDO 

c      WRITE(12,*)'in potential routine'
c      WRITE(12,*)'NMOLS =',NMOLS
c      DO I = 1, NMOLS
c         DO J = 1, NATOM(SPC(I))
c         WRITE(12,*)(R(I,J,IR),IR = 1,3)
c         END DO
c      END DO 
c      GOTO 444
c     THIS SECTION WILL BE ACCESSED FOR EVERY MONTE CARLO STEP
C     DUE TO THE FACT THAT ALL MOLECULES ARE MOVED AT ONCE
            
      DO I = 1, NMOLS - 1

      DO J = I + 1, NMOLS

            DO II = 1, NATOM(SPC(I)) 

                DO JJ = 1, NATOM(SPC(J)) 
    
                RXIIJJ = R(I,II,1) - R(J,JJ,1)
                RYIIJJ = R(I,II,2) - R(J,JJ,2)
                RZIIJJ = R(I,II,3) - R(J,JJ,3)
                RIIJJSQ = RXIIJJ**2 + RYIIJJ**2 + RZIIJJ**2
                RIIJJ = DSQRT(RIIJJSQ)

c             WRITE(12,*)'I =',I,' and J =',J
c             WRITE(12,*)'II =',II,' and JJ =',JJ
c             WRITE(12,*)'RIIJJ =',RIIJJ

C              This checks for water - water pairs
C	       IF(ANDX(SPC(I),SPC(J)).EQ.4)THEN

                IF(II.EQ.4)THEN
       
                PLRZ = J
                PLRZA = JJ
    
                RXIIJJ = R(I,II,1) - R(J,JJ,1)
                RYIIJJ = R(I,II,2) - R(J,JJ,2)
                RZIIJJ = R(I,II,3) - R(J,JJ,3)
                RIIJJSQ = RXIIJJ**2 + RYIIJJ**2 + RZIIJJ**2
                RIIJJ = DSQRT(RIIJJSQ)

c	        WRITE(12,*)'RIIJJ =',RIIJJ,' ANGSTROMS'

                EFUX = RXIIJJ/RIIJJ
                EFUY = RYIIJJ/RIIJJ
                EFUZ = RZIIJJ/RIIJJ
        
                EFLDMAG = dntmc_EFMAG(PLRZ,PLRZA,RIIJJSQ)
      
                EF(I,J,II,JJ,1) = EFLDMAG*EFUX
                EF(I,J,II,JJ,2) = EFLDMAG*EFUY
                EF(I,J,II,JJ,3) = EFLDMAG*EFUZ
        
                EFWX(I) = EFWX(I) + EF(I,J,II,JJ,1)
     +                    *ELAM(LNDX(I),LNDX(J))

                EFWY(I) = EFWY(I) + EF(I,J,II,JJ,2)
     +                    *ELAM(LNDX(I),LNDX(J))

                EFWZ(I) = EFWZ(I) + EF(I,J,II,JJ,3)
     +                    *ELAM(LNDX(I),LNDX(J))

c                WRITE(12,*)'EF(',I,',',J,',',II,',',JJ,',',1,') =',
c     +          EF(I,J,II,JJ,1)	
c                WRITE(12,*)'EF(',I,',',J,',',II,',',JJ,',',2,') =',
c     +          EF(I,J,II,JJ,2)
c                WRITE(12,*)'EF(',I,',',J,',',II,',',JJ,',',3,') =',
c     +          EF(I,J,II,JJ,3)

c                WRITE(12,*)'ELAM(',LNDX(I),',',LNDX(J),') =',
c     +          ELAM(LNDX(I),LNDX(J))

c               WRITE(12,*)'EFWX(',I,') =',EFWX(I)
c               WRITE(12,*)'EFWY(',I,') =',EFWY(I)
c               WRITE(12,*)'EFWZ(',I,') =',EFWZ(I)

                 IF((II.EQ.4).AND.(JJ.EQ.4))THEN

                   EF(J,I,JJ,II,1) = -EF(I,J,II,JJ,1)
                   EF(J,I,JJ,II,2) = -EF(I,J,II,JJ,2)
                   EF(J,I,JJ,II,3) = -EF(I,J,II,JJ,3)
           
                   EFWX(J) = EFWX(J) + EF(J,I,JJ,II,1)
     +                    *ELAM(LNDX(I),LNDX(J))

                   EFWY(J) = EFWY(J) + EF(J,I,JJ,II,2)
     +                    *ELAM(LNDX(I),LNDX(J))

                   EFWZ(J) = EFWZ(J) + EF(J,I,JJ,II,3)
     +                    *ELAM(LNDX(I),LNDX(J))


                   NUX(I,J) = RXIIJJ/RIIJJ
                   NUY(I,J) = RYIIJJ/RIIJJ
                   NUZ(I,J) = RZIIJJ/RIIJJ
c234567                   
C                   CHECKU=DSQRT(NUX(I,J)**2+NUY(I,J)**2
C     +             +NUZ(I,J)**2)

C                   WRITE(12,*)'II =',II,' and JJ =',JJ
C                   WRITE(12,*)'NUX(',I,',',J,') =',NUX(I,J)
C                   WRITE(12,*)'NUY(',I,',',J,') =',NUY(I,J)
C                   WRITE(12,*)'NUZ(',I,',',J,') =',NUZ(I,J)
C                   WRITE(12,*)'check for unity =',CHECKU

                   NUX(J,I) = -NUX(I,J)
                   NUY(J,I) = -NUY(I,J)
                   NUZ(J,I) = -NUZ(I,J)

C                   CHECKU=DSQRT(NUX(J,I)**2+NUY(J,I)**2
C     +             +NUZ(J,I)**2)

C                   WRITE(12,*)'JJ =',JJ,' and II =',II
C                   WRITE(12,*)'NUX(',J,',',I,') =',NUX(J,I)
C                   WRITE(12,*)'NUY(',J,',',I,') =',NUY(J,I)
C                   WRITE(12,*)'NUZ(',J,',',I,') =',NUZ(J,I)
C                   WRITE(12,*)'check for unity =',CHECKU

                   R3(I,J) = RIIJJ*RIIJJ*RIIJJ*6.748332695D+00
                   R3(J,I) = R3(I,J)
                   
C                   WRITE(12,*)'R3(',I,',',J,') =',R3(I,J)
C                   WRITE(12,*)'R3(',J,',',I,') =',R3(J,I)

c23456 	           
c                WRITE(12,*)'EF(',J,',',I,',',JJ,',',II,',',1,') =',
c     +          EF(J,I,JJ,II,1)	
c                WRITE(12,*)'EF(',J,',',I,',',JJ,',',II,',',2,') =',
c     +          EF(J,I,JJ,II,2)
c                WRITE(12,*)'EF(',J,',',I,',',JJ,',',II,',',3,') =',
c     +          EF(J,I,JJ,II,3)

c                WRITE(12,*)'ELAM(',LNDX(I),',',LNDX(J),') =',
c     +          ELAM(LNDX(I),LNDX(J))

c               WRITE(12,*)'EFWX(',J,') =',EFWX(J)
c               WRITE(12,*)'EFWY(',J,') =',EFWY(J)
c               WRITE(12,*)'EFWZ(',J,') =',EFWZ(J)
           
                 ENDIF

                ELSEIF(JJ.EQ.4)THEN
                
                PLRZ = I
                PLRZA = II
    
                RXIIJJ = R(J,JJ,1) - R(I,II,1)
                RYIIJJ = R(J,JJ,2) - R(I,II,2)
                RZIIJJ = R(J,JJ,3) - R(I,II,3)
                RIIJJSQ = RXIIJJ**2 + RYIIJJ**2 + RZIIJJ**2
                RIIJJ = DSQRT(RIIJJSQ)
        
C	        WRITE(6,*)'RIIJJ =',RIIJJ,' ANGSTROMS'

                EFUX = RXIIJJ/RIIJJ
                EFUY = RYIIJJ/RIIJJ
                EFUZ = RZIIJJ/RIIJJ
        
                EFLDMAG = dntmc_EFMAG(PLRZ,PLRZA,RIIJJSQ)
      
                EF(J,I,JJ,II,1) = EFLDMAG*EFUX
                EF(J,I,JJ,II,2) = EFLDMAG*EFUY
                EF(J,I,JJ,II,3) = EFLDMAG*EFUZ
        
                EFWX(J) = EFWX(J) + EF(J,I,JJ,II,1)
     +                    *ELAM(LNDX(I),LNDX(J))

                EFWY(J) = EFWY(J) + EF(J,I,JJ,II,2)
     +                    *ELAM(LNDX(I),LNDX(J))

                EFWZ(J) = EFWZ(J) + EF(J,I,JJ,II,3)
     +                    *ELAM(LNDX(I),LNDX(J))

        
c                WRITE(12,*)'EF(',J,',',I,',',JJ,',',II,',',1,') =',
c     +          EF(J,I,JJ,II,1)	
c                WRITE(12,*)'EF(',J,',',I,',',JJ,',',II,',',2,') =',
c     +          EF(J,I,JJ,II,2)
c                WRITE(12,*)'EF(',J,',',I,',',JJ,',',II,',',3,') =',
c     +          EF(J,I,JJ,II,3)

c                WRITE(12,*)'ELAM(',LNDX(I),',',LNDX(J),') =',
c     +          ELAM(LNDX(I),LNDX(J))

c                WRITE(12,*)'EFWX(',J,') =',EFWX(J)
c                WRITE(12,*)'EFWY(',J,') =',EFWY(J)
c                WRITE(12,*)'EFWZ(',J,') =',EFWZ(J)
     
                ENDIF
C	       ENDIF
     

                                                          
               IF(II.EQ.1.AND.JJ.EQ.1)THEN

                 Vh2o(II,JJ) = dntmc_Vljones(I,J,II,JJ,RIIJJSQ)
c                WRITE(12,*)'Vh2oOO =',Vh2o(II,JJ),' kcal/mole'              
                                               
               ELSEIF(II.EQ.4.AND.JJ.EQ.4)THEN
       
                 Vh2o(II,JJ) = dntmc_Vcoul(I,J,II,JJ,RIIJJ)
c                WRITE(12,*)'Vh2o =',Vh2o(II,JJ)
                 
               ELSEIF(II.EQ.1.OR.JJ.EQ.1)THEN

                 Vh2o(II,JJ) = 0.D+00               
                 
               ELSE
               
                 Vh2o(II,JJ) = dntmc_Vcoul(I,J,II,JJ,RIIJJ)
c                WRITE(12,*)'Vh2o =',dntmc_Vcoul(I,J,II,JJ,RIIJJ)              

               ENDIF


               V(I,J,II,JJ) = Vh2o(II,JJ)      

c          WRITE(12,*)'LNDX(',I,') =',LNDX(I)
c          WRITE(12,*)'LNDX(',J,') =',LNDX(J)

c          WRITE(12,*)'V(',I,',',J,',',II,',',JJ,') =',V(I,J,II,JJ)
c234567                	       
               UU = UU + V(I,J,II,JJ)*
     +              LAM(AB,LNDX(I),LNDX(J),LNDXX(SPC(I),SPC(J)))

    
               END DO  ! JJ

            END DO     ! II  
 
         END DO        ! J

      END DO           ! I

c      WRITE(12,*)'UU and DU before induced dipoles are computed'
c      WRITE(12,*)'UU =',UU

C      WRITE(12,*)'Done with pair potentials and point charge'
c      WRITE(12,*)'electric field calculation...'
C      WRITE(12,*)'...'

      
C This section calculates the self-consistent induced dipole moments
C for each molecule and once converged evaluates Upol

C This first loop is over the number of iterations required to
C achieve self-consistency

      POL(2) = 9.74459258D+00
      MUTOL = 1.D-6

      DO I = 1, NMOLS
         COUNT(I) = 0.D+00
         MU(1,I) = 0.D+00
         MUX(1,I) = 0.D+00
         MUY(1,I) = 0.D+00
         MUZ(1,I) = 0.D+00
         EFX(I) = 0.D+00
         EFY(I) = 0.D+00
         EFZ(I) = 0.D+00
      ENDDO

      SUM = 0.D+00
c      WRITE(12,*)'Beginning self-consistent dipole calculation'
c      WRITE(12,*)'first loop...'

      DO SCF = 2, NSCF

C        WRITE(12,*)'SCF =',SCF

        DO I = 1, NMOLS

          DO J = 1, I - 1

          DOTX = NUX(I,J)*MUX(SCF-1,J)
          DOTY = NUY(I,J)*MUY(SCF-1,J)
          DOTZ = NUZ(I,J)*MUZ(SCF-1,J)

          DOT = DOTX + DOTY + DOTZ

          EFIX(I)=(3.D+00*DOT*NUX(I,J)-MUX(SCF-1,J))/R3(I,J)
          EFIY(I)=(3.D+00*DOT*NUY(I,J)-MUY(SCF-1,J))/R3(I,J)
          EFIZ(I)=(3.D+00*DOT*NUZ(I,J)-MUZ(SCF-1,J))/R3(I,J)

          EFIX(I) = EFIX(I)*ELAM(LNDX(I),LNDX(J))
          EFIY(I) = EFIY(I)*ELAM(LNDX(I),LNDX(J))
          EFIZ(I) = EFIZ(I)*ELAM(LNDX(I),LNDX(J))

          EFX(I) = EFX(I) + EFIX(I)
          EFY(I) = EFY(I) + EFIY(I)
          EFZ(I) = EFZ(I) + EFIZ(I)

          ENDDO

          DO J = I + 1, NMOLS

          DOTX = NUX(I,J)*MUX(SCF-1,J)
          DOTY = NUY(I,J)*MUY(SCF-1,J)
          DOTZ = NUZ(I,J)*MUZ(SCF-1,J)

          DOT = DOTX + DOTY + DOTZ

          EFIX(I)=(3.D+00*DOT*NUX(I,J)-MUX(SCF-1,J))/R3(I,J)
          EFIY(I)=(3.D+00*DOT*NUY(I,J)-MUY(SCF-1,J))/R3(I,J)
          EFIZ(I)=(3.D+00*DOT*NUZ(I,J)-MUZ(SCF-1,J))/R3(I,J)

          EFIX(I) = EFIX(I)*ELAM(LNDX(I),LNDX(J))
          EFIY(I) = EFIY(I)*ELAM(LNDX(I),LNDX(J))
          EFIZ(I) = EFIZ(I)*ELAM(LNDX(I),LNDX(J))

          EFX(I) = EFX(I) + EFIX(I)
          EFY(I) = EFY(I) + EFIY(I)
          EFZ(I) = EFZ(I) + EFIZ(I)
          
          ENDDO

        ENDDO  ! I loop

C        WRITE(12,*)'in 2nd I loop estimating induced moments...'

c23456789

        DO I = 1, NMOLS

          MUX(SCF,I)=POL(SPC(I))*(EFWX(I) + EFX(I))
          MUY(SCF,I)=POL(SPC(I))*(EFWY(I) + EFY(I))
          MUZ(SCF,I)=POL(SPC(I))*(EFWZ(I) + EFZ(I))
          
          MU(SCF,I)=DSQRT(MUX(SCF,I)**2+MUY(SCF,I)**2+MUZ(SCF,I)**2)

          DMU=DABS(MU(SCF,I)-MU(SCF-1,I))

          IF(DMU.LT.MUTOL)THEN

            MX(I)= MUX(SCF,I)
            MY(I)= MUY(SCF,I)
            MZ(I)= MUZ(SCF,I)

            X=DSQRT(MX(I)**2+MY(I)**2+MZ(I)**2)

c            WRITE(12,*)'hey dipole',I,' converged'
c            WRITE(12,*)'after ',SCF,' steps'
c            WRITE(12,*)'************************'
c            WRITE(12,*)'MX(',I,') =',MX(I)
c            WRITE(12,*)'MY(',I,') =',MY(I)
c            WRITE(12,*)'MZ(',I,') =',MZ(I)
c            WRITE(12,*)'|MU| =',X
c            WRITE(12,*)'************************'

            COUNT(I) = 1.D+00

C            WRITE(12,*)'COUNT(',I,') =',COUNT(I)            

          ENDIF

        ENDDO  ! I
 
        DO J = 1, NMOLS
           SUM = SUM + COUNT(J)
        ENDDO

          
        IF(SUM.EQ.DBLE(NMOLS))THEN

C          WRITE(12,*)'all induced dipoles have converged!'
C          WRITE(12,*)'SUM =',SUM
          
          GOTO 817

        ENDIF 
           
        SUM = 0.D+00

        DO J = 1, NMOLS
           EFX(J) = 0.D+00
           EFY(J) = 0.D+00
           EFZ(J) = 0.D+00
           COUNT(J) = 0.D+00
        ENDDO

C        WRITE(12,*)'finished with lower ith loop'
C        WRITE(12,*)'going back for another SCF iteration'

        IF(SCF.GE.NSCF)THEN
    
        WRITE(12,*)'...there is a problem...'
        WRITE(12,*)'the induced dipoles did not converge'
        WRITE(12,*)'...stopping...sorry'

        call errquit('dntmc_psum: did not converge', 1, CALC_ERR)

        ENDIF

      ENDDO   ! SCF loop

 817  CONTINUE

C******************************************************************
C This loop sums the dot products of the self consistent induced
C dipole moments with the total electric field of the charges
C to give the polarization interaction

      DO I = 1, NMOLS
         DOTX=MX(I)*EFWX(I)
         DOTY=MY(I)*EFWY(I)
         DOTZ=MZ(I)*EFWZ(I)

         DOT = DOTX + DOTY + DOTZ
         UDOT1 = - 0.5D+00*DOT*627.5095D+00
         Upol1 = Upol1 + UDOT1

c         WRITE(12,*)'UDOT1(',I,') =',UDOT1
c         WRITE(12,*)'Upol1(',I,') =',Upol1
      ENDDO

c      WRITE(12,*)'before these energies are added...'
c      WRITE(12,*)'UU =',UU,' kcal/mole'
c      WRITE(12,*)'Upol1 =',Upol1,' kcal/mole'
c      WRITE(12,*)'AB =',AB

c      WRITE(12,*)'these are the total dipole moments...'

c      DO I = 1, NMOLS
c	RXIIJJ = R(I,1,1) - R(I,4,1)
c	RYIIJJ = R(I,1,2) - R(I,4,2)
c	RZIIJJ = R(I,1,3) - R(I,4,3)
c	RIIJJSQ = RXIIJJ**2 + RYIIJJ**2 + RZIIJJ**2
c	RIIJJ = DSQRT(RIIJJSQ)

c        MNUX = -0.7279259771D+00*RXIIJJ/RIIJJ
c        MNUY = -0.7279259771D+00*RYIIJJ/RIIJJ
c        MNUZ = -0.7279259771D+00*RZIIJJ/RIIJJ

c        MTX = MNUX + MX(I) 
c        MTY = MNUY + MY(I) 
c        MTZ = MNUZ + MZ(I) 

c        MTMAG = DSQRT(MTX**2 + MTY**2 + MTZ**2)

c        MTMAG = MTMAG*2.541467207D+00

c        WRITE(12,*)'MUtot(',I,') =',MTMAG,' Debye'
c      ENDDO

      UU = UU + Upol1

c 444  CONTINUE
c      UU = 0.D+00
c      WRITE(6,*)'leaving psum'

      RETURN
      END

C***********************************************************************
      double precision FUNCTION dntmc_RAN1(IDUM)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
c     N.R. FUNCTION ran2(idum) eventhough it is called RAN1
c
      IMPLICIT none

      integer idum   ! input

      integer im1, im2, imm1, ia1, ia2, iq1, iq2, ir1, ir2, ntab, ndiv
      double precision rnmx
      double precision AM,EPS

      PARAMETER (IM1 = 2147483563)
      PARAMETER (IM2 = 2147483399)
      PARAMETER (AM = 1.D+00/IM1)
      PARAMETER (IMM1 = IM1 - 1)
      PARAMETER (IA1 = 40014)
      PARAMETER (IA2 = 40692)
      PARAMETER (IQ1 = 53668)
      PARAMETER (IQ2 = 52774)
      PARAMETER (IR1 = 12211)
      PARAMETER (IR2 = 3791)
      PARAMETER (NTAB = 32)
      PARAMETER (NDIV = 1 + IMM1/NTAB)
      PARAMETER (EPS = 1.2D-07)
      PARAMETER (RNMX = 1.D+00 - EPS)

      INTEGER IV(NTAB)
      integer iy, idum2, j, k

      SAVE IV,IY,IDUM2

      DATA IDUM2/123456789/, IV/NTAB*0/, IY/0/

      IF (IDUM.LE.0) THEN

        IDUM = MAX(-IDUM,1)
        IDUM2 = IDUM

        DO J = NTAB + 8, 1, -1
           K = IDUM/IQ1
           IDUM = IA1*(IDUM - K*IQ1) - K*IR1

           IF (IDUM.LT.0) IDUM = IDUM + IM1
           IF (J.LE.NTAB) IV(J) = IDUM
        END DO

        IY = IV(1)

      ENDIF

      K = IDUM/IQ1

      IDUM = IA1*(IDUM - K*IQ1) - K*IR1
   
      IF (IDUM.LT.0) IDUM = IDUM + IM1

      K = IDUM2/IQ2

      IDUM2 = IA2*(IDUM2 - K*IQ2) - K*IR2

      IF (IDUM2.LT.0) IDUM2 = IDUM2 + IM2

      J = 1 + IY/NDIV

      IY = IV(J) - IDUM2

      IV(J) = IDUM
    
      IF (IY.LT.1) IY = IY + IMM1

      dntmc_RAN1 = MIN(AM*IY,RNMX)

      RETURN
      END

C***********************************************************************
      double precision FUNCTION dntmc_RAN3(IDUM)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
c     N.R. FUNCTION ran3(idum)
c
      implicit none

      integer idum  ! input

      integer mbig, mseed, ma
      integer iff, mj, mk, mz, inext, inextp
      integer i, ii, k
      double precision fac
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      dimension ma(55)
      data iff /0/
      save iff, mj, ma, mk, inext, inextp
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      dntmc_ran3=mj*fac
      return
      end

C***********************************************************************
      SUBROUTINE dntmc_RAN(KK)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none

      integer kk ! input

      double precision seed
      double precision rnnum

      DOUBLE PRECISION A,B

      COMMON/SEEDS/SEED(1)
      COMMON/RNNUMS/RNNUM(1)

      A = (7.0D+00**5)*SEED(KK)
      B = 2.0D+00**31 - 1.D+00
      SEED(KK) = DMOD(A,B)
      RNNUM(KK) = SEED(KK) / 2.0D+00**31

      RETURN
      END

C***********************************************************************
      double precision FUNCTION dntmc_Vcoul(I,J,II,JJ,RIIJJ)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"

      integer i, j, ii, jj ! input
      double precision riijj ! input

      double precision qsqrd, sgmsq, epsln
      INTEGER SPC
      double precision e, k, esqrd, kb, pi
      integer l

      COMMON/INDEXX/QSQRD(MXSP,MXSP,MXATOM,MXATOM),
     +SGMSQ(MXSP,MXSP,MXATOM,MXATOM),
     +EPSLN(MXSP,MXSP,MXATOM,MXATOM)
      COMMON/SPECIES/SPC(MXMOLS)     
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/MSTP/L
                      
      dntmc_Vcoul = (k*QSQRD(SPC(I),SPC(J),II,JJ)*esqrd) / RIIJJ
      
      RETURN
      END 

C***********************************************************************
      double precision FUNCTION dntmc_Vljones(I,J,II,JJ,RIIJJSQ)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"

      integer i, j, ii, jj  ! input
      double precision riijjsq ! input

      double precision qsqrd, sgmsq, epsln
      INTEGER SPC
      integer l

      double precision sr2, sr6, sr12

      COMMON/INDEXX/QSQRD(MXSP,MXSP,MXATOM,MXATOM),
     +SGMSQ(MXSP,MXSP,MXATOM,MXATOM),
     +EPSLN(MXSP,MXSP,MXATOM,MXATOM)
      COMMON/SPECIES/SPC(MXMOLS)
      COMMON/MSTP/L
                 
      SR2 = SGMSQ(SPC(I),SPC(J),II,JJ) / RIIJJSQ

      SR6 = SR2*SR2*SR2

      SR12 = SR6*SR6
      
c      PRINT*,'SR2 =',SR2
c      PRINT*,'SR6 =',SR6
c      PRINT*,'SR12 =',SR12
           
      dntmc_Vljones = 4.D+00*EPSLN(SPC(I),SPC(J),II,JJ)*(SR12 - SR6)
     
      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_AVERAGE(L,TEMP,UU,UAVR,Cv,RSIMAVR)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none

      integer l     ! input
      double precision temp, uu ! input
      double precision uavr, cv, rsimavr ! output

      double precision e, k, esqrd, kb, pi
      double precision usum, usqrdsum
      double precision rmin, rsim, rmax

      double precision rsimsum, usqrd, tsqrd

      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/AVR/USUM,USQRDSUM      
      COMMON/RADIUS/RMIN,RSIM,RMAX      

      SAVE RSIMSUM

      IF (L.EQ.1) THEN
         USUM = 0.D+00
         RSIMSUM = 0.D+00
         USQRDSUM = 0.D+00
      END IF
      
c     WRITE(12,*)'RSIM =',RSIM,'at step L =',L
      
      USUM = USUM + UU

      RSIMSUM = RSIMSUM + RSIM

      USQRDSUM = USQRDSUM + UU*UU

      UAVR = USUM/DBLE(L)

      RSIMAVR = RSIMSUM/DBLE(L)
         
      USQRD = USQRDSUM/DBLE(L)
      TSQRD = TEMP*TEMP

      Cv = (USQRD - (UAVR*UAVR))/(kb*TSQRD)

      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_ROOTMEANSQUARE(L,NMOLS,NSPC)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"
      integer l, nmols, nspc ! input - nspc is not used

      double precision rmsdsum, rmsdavr
      double precision r, rcm, rs
      double precision ri
      double precision rmsdsumsp, rmsdavrsp
      INTEGER SPC,USPC

      double precision rmsd, rmsdx, rmsdy, rmsdz
      integer i

      COMMON/RMSDCOORD/RMSDSUM(MXMOLS),RMSDAVR(MXMOLS)
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/COORDI/RI(MXMOLS,MXATOM,3)
      COMMON/RMSDCOORDSP/RMSDSUMSP(MXSP),RMSDAVRSP(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)
      COMMON/USPECIES/USPC(MXSP)
      
      DIMENSION RMSD(MXMOLS)
      
      DO I = 1, NMOLS
         RMSDX = ( R(I,1,1) - RI(I,1,1) )**2
         RMSDY = ( R(I,1,2) - RI(I,1,2) )**2
         RMSDZ = ( R(I,1,3) - RI(I,1,3) )**2
         RMSD(I) = RMSDX + RMSDY + RMSDZ
      
         RMSDSUM(I) = RMSDSUM(I) + RMSD(I)
      
         RMSDAVR(I) = RMSDSUM(I)/DBLE(L)
      
         IF(SPC(I).EQ.2)THEN
            RMSDSUMSP(SPC(I)) = RMSDSUMSP(SPC(I)) + RMSD(I) 
            RMSDAVRSP(SPC(I)) = RMSDSUMSP(SPC(I))/DBLE(L)
         ELSEIF(SPC(I).EQ.5)THEN
            RMSDSUMSP(SPC(I)) = RMSDSUMSP(SPC(I)) + RMSD(I) 
            RMSDAVRSP(SPC(I)) = RMSDSUMSP(SPC(I))/DBLE(L)
         ELSE
            RMSDSUMSP(SPC(I)) = RMSDSUMSP(SPC(I)) + RMSD(I) 
            RMSDAVRSP(SPC(I)) = RMSDSUMSP(SPC(I))/DBLE(L)
         ENDIF
      END DO      
            
      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_CLUSTERCMASS(NMOLS) 
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"

      integer nmols

      double precision M,MTOT
      double precision r, rcm, rs
      double precision rclcm
      INTEGER SPC

      double precision MCLTOT
      integer i

      COMMON/MASSES/M(MXSP,MXATOM),MTOT(MXSP)
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/COORD2/RCLCM(3)
      COMMON/SPECIES/SPC(MXMOLS)


      MCLTOT = 0.D+00
      RCLCM(1) = 0.D+00
      RCLCM(2) = 0.D+00
      RCLCM(3) = 0.D+00

      DO I = 1, NMOLS
         MCLTOT = MCLTOT + MTOT(SPC(I))

         CALL dntmc_CMASS(I)

         RCLCM(1) = MTOT(SPC(I))*RCM(I,1) + RCLCM(1)
         RCLCM(2) = MTOT(SPC(I))*RCM(I,2) + RCLCM(2)
         RCLCM(3) = MTOT(SPC(I))*RCM(I,3) + RCLCM(3)
      END DO

      RCLCM(1) = RCLCM(1)/MCLTOT
      RCLCM(2) = RCLCM(2)/MCLTOT
      RCLCM(3) = RCLCM(3)/MCLTOT


c      DO I = 1, 3
c         WRITE(12,*)'RCLCM(',I,') =',RCLCM(I)
c      END DO
  
      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_CMASS(I)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"

      integer i  ! input

      double precision M,MTOT
      double precision r, rcm, rs
      double precision e, k, esqrd, kb, pi
      integer natom
      INTEGER SPC

      integer j

      COMMON/MASSES/M(MXSP,MXATOM),MTOT(MXSP)
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)

      RCM(I,1) = 0.D+00
      RCM(I,2) = 0.D+00
      RCM(I,3) = 0.D+00

c      WRITE(12,*)'in cm routine'
c      WRITE(12,*)'MTOT(',SPC(I),') =',MTOT(SPC(I))

      DO J = 1, NATOM(SPC(I))
c          WRITE(12,*)'M(',SPC(I),',',J,') =',M(SPC(I),J)
          RCM(I,1)  = M(SPC(I),J)*R(I,J,1) + RCM(I,1)
          RCM(I,2)  = M(SPC(I),J)*R(I,J,2) + RCM(I,2)
          RCM(I,3)  = M(SPC(I),J)*R(I,J,3) + RCM(I,3)
      END DO

      RCM(I,1) = RCM(I,1) / MTOT(SPC(I))
      RCM(I,2) = RCM(I,2) / MTOT(SPC(I))
      RCM(I,3) = RCM(I,3) / MTOT(SPC(I))


c      DO J = 1, 3
c       WRITE(12,*)'RCM(',I,',',J,') =',RCM(I,J)
c      END DO

      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_NUMBERHIST(L)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none

      integer l  ! input

      double precision rmin, rsim, rmax
      integer number

      integer nob, i
      double precision INCR

      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/NUMBER1/NUMBER(100)    

      NOB = 100

      INCR = (RSIM - RMIN)/(RMAX - RMIN)

      IF(L.EQ.1)THEN
         DO I = 1, NOB
            NUMBER(I) = 0
         END DO
      endif

      I = NINT(INCR*100.D+00)

c      WRITE(12,*)'RSIM =',RSIM
c      WRITE(12,*)'INCR =',INCR
c      WRITE(12,*)'I =',I

        NUMBER(I) = NUMBER(I) + 1

c      WRITE(12,*)'NUMBER(',I,') =',NUMBER(I)

      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_NUMBERAVERAGE(L,temp,rweight)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none

      integer l   ! input

      integer number
      double precision NUMBERAVR,NUMBERSUM
      double precision rb, testme
      integer nob, i
      double precision e,k,esqrd,kb,PI    
      double precision temp
      logical rweight 

      COMMON/NUMBER1/NUMBER(100)
      COMMON/NUMBER2/NUMBERSUM(100),NUMBERAVR(100)
      COMMON/DIST1/RB(100),TESTME(100,2000)
      COMMON/CONSTANTS/e,k,esqrd,kb,PI

      NOB = 100

      IF(L.EQ.1)THEN
        DO I = 1, NOB
           NUMBERSUM(I) = 0.D+00
           NUMBERAVR(I) = 0.D+00
        END DO
      ENDIF

      if (rweight) then
        DO I = 1, NOB
          NUMBERSUM(I) = DBLE(NUMBER(I))/RB(I)**2.D+00
          if (i .gt. 1) then
            NUMBERSUM(I) = NUMBERSUM(I) + NUMBERSUM(I-1)
          endif
        enddo
        DO I = 1, NOB
          NUMBERSUM(I) = NUMBERSUM(I)/NUMBERSUM(NOB)
          NUMBERSUM(I) = NUMBERSUM(I)*4.D+00*PI*RB(I)**2.D+00
          NUMBERSUM(I) = NUMBERSUM(I)*DEXP(-1.7D-05*4.D+00/3.D+00*PI
     &    *(RB(I)**3.D+00 - RB(NOB)**3.D+00)/kb/TEMP)
          NUMBERAVR(I) = NUMBERSUM(I)
        ENDDO
      else   
      DO I = 1, NOB
C  Don't believe this second sum is needed
C         NUMBERSUM(I) = NUMBERSUM(I) + DBLE(NUMBER(I))
         NUMBERSUM(I) = DBLE(NUMBER(I))
C  End Corrected NUMERSUM
         NUMBERAVR(I) = NUMBERSUM(I)/DBLE(L)
      END DO
      endif

      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_ANGLEAVR(L)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none

      integer L   ! input

      double precision phi, phisum, phiavr

      integer i

      COMMON/ANGLES/PHI(100),PHISUM(100),PHIAVR(100)

c     IF(L.EQ.1)THEN

      DO I = 1,100
         PHISUM(I) = 0.D+00
         PHIAVR(I) = 0.D+00
      END DO

c     ENDIF

      DO I = 1,100
         PHISUM(I) = PHISUM(I) + PHI(I)
         PHIAVR(I) = PHISUM(I)/DBLE(L)
      END DO

      RETURN
      END
      
C***********************************************************************
      double precision FUNCTION dntmc_EFMAG(PLRZ,PLRZA,RIIJJSQ)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"

      integer plrz, plrza       ! input
      double precision riijjsq  ! input

      double precision qse
      INTEGER SPC
      double precision e, k, esqrd, kb, pi
      integer l

      integer i, ii

      COMMON/INDEX/QSE(MXSP,MXATOM,3)
      COMMON/SPECIES/SPC(MXMOLS)     
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/MSTP/L

      I = PLRZ
      II = PLRZA

C      WRITE(6,*)'in dntmc_efmag function call'
C      WRITE(6,*)'PLRZA =',II
C      WRITE(6,*)'Q(',II,') =',QSE(SPC(I),II,1)

      dntmc_EFMAG= (k*e*QSE(SPC(I),II,1))/(RIIJJSQ)
      dntmc_EFMAG= DBLE(dntmc_EFMAG/(7.401250726D+21))
      
c      WRITE(12,*)'Q =',QSE(SPC(I),II,1)
c      WRITE(12,*)'dntmc_EFMAG(',I,',',II,') =',dntmc_EFMAG,' E-FIELD'
      
      RETURN
      END         

C***********************************************************************
      SUBROUTINE dntmc_CONFIG(IMLW,NODE,RMX,rinit)
C
C Make sure molecules don't overlap
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"
      integer imlw, node    ! input
      double precision rmx  ! input

      double precision M,MTOT
      double precision r, rcm, rs
      double precision rc
      double precision pi, d
      integer natom
      double precision displ, angle, rmxsq
      INTEGER SPC

      integer USPC, mcsteps, i, j, k, ll
      double precision rinit, rnode, pct, pctt
      LOGICAL ACCP

      COMMON/MASSES/M(MXSP,MXATOM),MTOT(MXSP)
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/CORD/RC(MXMOLS,MXATOM,3)
      COMMON/CNSTNTS/PI,D
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/TRANSROT/DISPL,ANGLE,RMXSQ
      COMMON/SPECIES/SPC(MXMOLS)

      DIMENSION USPC(1),RINIT(MXMOLS,MXATOM,3)


      PI = 3.141592653589793D+00
      
      D = 2.5D+00  !A = minimal distance between atoms

      USPC(1) = 2

      WRITE(12,*)'INITIAL CONFIG DATA'
      WRITE(12,*)'IMLW =',IMLW
      WRITE(12,*)'NODE =',NODE
      WRITE(12,*)'RMX =',RMX

      DO I = 1, IMLW
         WRITE(12,*)'SPC(',I,') =',SPC(I)
      ENDDO

      WRITE(12,*)'NATOM(',USPC(1),') =',NATOM(USPC(1))
            
      MCSTEPS = 100   !MCSTEPS

      RNODE = DBLE(NODE*NODE)/938.D+00 

      PCT = RNODE*1.D+00 + 25.D+00   !PCT  % of RMAX
      PCTT = RNODE*1.D+00 + 30.D+00  !PCTT  % of RMAX^2

      ANGLE = 1.0D+00*PI
 
      RMXSQ = RMX*RMX
      
      RMXSQ = (DBLE(PCTT)/100.D+00)*RMXSQ
                                      
      DISPL = (DBLE(PCT)/100.D+00)*RMX
     
C*********************************************************************
C     COORDINATES FOR DANG-CHANG WATER MODEL WITH C2v in Z-DIRECTION
c this will need to go away since the initial coordinates should
c come from the nwchem input, it should be passed in
c     RINIT(1,1,1) = 0.D+00
c     RINIT(1,1,2) = 0.D+00
c     RINIT(1,1,3) = 0.D+00
c     RINIT(1,2,1) = 0.D+00
c     RINIT(1,2,2) = 0.75695D+00
c     RINIT(1,2,3) = 0.58588D+00
c     RINIT(1,3,1) = 0.D+00
c     RINIT(1,3,2) = -0.75695D+00
c     RINIT(1,3,3) = 0.58588D+00
c     RINIT(1,4,1) = 0.D+00
c     RINIT(1,4,2) = 0.D+00
c     RINIT(1,4,3) = 0.2150D+00
C*********************************************************************
      
c      WRITE(12,*)'coords before moving initial molecules:'

c needs to be rinit(spc(i),j,1..3)?
      DO I = 1, IMLW
         DO J = 1, NATOM(SPC(I))
            RC(I,J,1) = RINIT(1,J,1)
            RC(I,J,2) = RINIT(1,J,2)
            RC(I,J,3) = RINIT(1,J,3)
c            WRITE(12,*)(RC(I,J,IR),IR=1,3)
         END DO
      END DO
      
      
      DO K = 1, IMLW

      WRITE(12,*)'MOLECULE =',K
       
      DO LL = 1 , MCSTEPS
         ACCP = .FALSE.
         
         CALL dntmc_MOVEM(K,IMLW,ACCP)
         
         IF(ACCP) goto 77
      END DO
      
  77  CONTINUE
  
      WRITE(12,*)'MOLECULE =',K,' WAS MOVED ',LL,' TIMES'
      
      END DO ! K

c      WRITE(12,*)'****************************************'
c      WRITE(12,*)'coordinates from config subroutine...'

      DO I = 1, IMLW
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = RC(I,J,1)
            R(I,J,2) = RC(I,J,2)
            R(I,J,3) = RC(I,J,3)
c	    WRITE(12,*)'RC(',I,',',J,',',1,') =',RC(I,J,1)
c	    WRITE(12,*)'RC(',I,',',J,',',2,') =',RC(I,J,2)
c	    WRITE(12,*)'RC(',I,',',J,',',3,') =',RC(I,J,3)
c	    WRITE(12,*)'R(',I,',',J,',',1,') =',R(I,J,1)
c	    WRITE(12,*)'R(',I,',',J,',',2,') =',R(I,J,2)
c	    WRITE(12,*)'R(',I,',',J,',',3,') =',R(I,J,3)
         END DO
      END DO

c      WRITE(12,*)'****************************************'

      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_MOVEM(K,IMLW,ACCP)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"
      integer k, imlw   ! input
      LOGICAL ACCP      ! output

      double precision r, rcm, rs
      double precision rc
      integer natom
      double precision displ, angle, rmxsq
      INTEGER SPC
      double precision seed
      double precision rnnum
      double precision dlta, rsave

      double precision comp2, comp3, comp2r, comp3r
      double precision beta, sdbeta, cdbeta, rx, ry, rz, dist
      integer i, j, ia1, ia2, ia3

      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/CORD/RC(MXMOLS,MXATOM,3)
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/TRANSROT/DISPL,ANGLE,RMXSQ
      COMMON/SPECIES/SPC(MXMOLS)      
c      COMMON/SEED/IDUM 
      COMMON/SEEDS/SEED(1)
      COMMON/RNNUMS/RNNUM(1)

      DIMENSION DLTA(3),RSAVE(MXATOM,3)

      COMP2 = 0.D+00
      COMP3 = 0.D+00
      COMP2R = 0.D+00
      COMP3R = 0.D+00

C     SAVE THE COORDINATES OF THE MOLECULE K
  
      DO I = 1, NATOM(SPC(K)) 
           DO J = 1, 3
              RSAVE(I,J) = RC(K,I,J)
              WRITE(12,*)'RSAVE(',I,',',J,') =',RSAVE(I,J)
           END DO
      END DO

C     PICK DISPLACEMENTS ALONG X,Y,Z
      
      CALL dntmc_RAN(1)
      DLTA(1) = DISPL*(2.D+00*RNNUM(1) - 1.D+00)
      CALL dntmc_RAN(1)
      DLTA(2) = DISPL*(2.D+00*RNNUM(1) - 1.D+00)
      CALL dntmc_RAN(1)
      DLTA(3) = DISPL*(2.D+00*RNNUM(1) - 1.D+00)


      DO I = 1, 3
         WRITE(12,*)'DLTA(',I,') =',DLTA(I)
      ENDDO
      
C     TRANSLATE COORDINATES OF CHOSEN MOLECULE
c      WRITE(12,*)'coords after translation'

      DO I = 1, NATOM(SPC(K))
         DO J = 1, 3
            RC(K,I,J) = RC(K,I,J) + DLTA(J)
            WRITE(12,*)'RC(',K,',',I,',',J,') =',RC(K,I,J)
         END DO
      END DO

C     PICK AN ANGLE

      CALL dntmc_RAN(1)
      BETA = 2.D+00*ANGLE*(RNNUM(1) - 0.5D+00)

      WRITE(12,*)'ROTATION ANGLE =',BETA,' radians'

      SDBETA = DSIN(BETA)
      CDBETA = DCOS(BETA)

C     PICK AN AXIS ABOUT WHICH TO ROTATE

      CALL dntmc_RAN(1)
      IA1 = INT(3.D+00*RNNUM(1) + 1.D+00)
      IA2 = IA1 + 1
      IA3 = IA1 + 2
      IF(IA2.GT.3) IA2 = IA2 - 3
      IF(IA3.GT.3) IA3 = IA3 - 3

      WRITE(12,*)'ROT.AXIS(xyz->123) =',IA1,IA2,IA3

C     ROTATE THE CHOSEN MOLECULE ABOUT CENTER OF MASS
      
      DO J = 1, NATOM(SPC(K))
            R(K,J,1) = RC(K,J,1)
            R(K,J,2) = RC(K,J,2)
            R(K,J,3) = RC(K,J,3)
      END DO
      
      CALL dntmc_CMASS(K)

      DO J = 1, NATOM(SPC(K))
            RC(K,J,1) = R(K,J,1)
            RC(K,J,2) = R(K,J,2)
            RC(K,J,3) = R(K,J,3)
      END DO

      DO J = 1, NATOM(SPC(K))
         COMP2 = (RC(K,J,IA2) - RCM(K,IA2))
         COMP3 = (RC(K,J,IA3) - RCM(K,IA3))

         COMP2R=CDBETA*COMP2+SDBETA*COMP3
         COMP3R=-SDBETA*COMP2+CDBETA*COMP3

         RC(K,J,IA2)= RCM(K,IA2) + COMP2R
         RC(K,J,IA3)= RCM(K,IA3) + COMP3R

         RX = RC(K,J,IA1)
         RY = RC(K,J,IA2)
         RZ = RC(K,J,IA3)

         DIST = RX*RX + RY*RY + RZ*RZ         

         IF(DIST.GT.RMXSQ) GOTO 11
      END DO

      DO J = 1, NATOM(SPC(K))
         WRITE(12,*)'RC(',K,',',J,',',1,') =',RC(K,J,1)
         WRITE(12,*)'RC(',K,',',J,',',2,') =',RC(K,J,2)
         WRITE(12,*)'RC(',K,',',J,',',3,') =',RC(K,J,3)
      END DO

      CALL dntmc_STRUCT2(K,IMLW,ACCP)            
  
c      IF (ACCP) GOTO 12
      IF (ACCP) then
      GOTO 12
      
  11  CONTINUE
  
      ACCP = .FALSE.

      write(12,*) 'inside of 11 continue'
      
      DO I = 1, NATOM(SPC(K))
           DO J = 1, 3
               RC(K,I,J) = RSAVE(I,J)
           END DO
      END DO

      endif
      
  12  CONTINUE

      WRITE(12,*)'after structure routine'
      DO J = 1, NATOM(SPC(K))
         WRITE(12,*)'RC(',K,',',J,',',1,') =',RC(K,J,1)
         WRITE(12,*)'RC(',K,',',J,',',2,') =',RC(K,J,2)
         WRITE(12,*)'RC(',K,',',J,',',3,') =',RC(K,J,3)
      END DO

      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_STRUCT2(K,IMLW,ACCP)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"
      integer k, imlw    ! input
      LOGICAL ACCP       ! output

      double precision r, rcm, rs
      double precision rc
      double precision pi, d
      integer natom
      INTEGER SPC

      integer i, ii, jj
      double precision rxiijj, ryiijj, rziijj, riijjsq
      double precision riijj

      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/CORD/RC(MXMOLS,MXATOM,3)
      COMMON/CNSTNTS/PI,D
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)

      WRITE(12,*)'d is ', d

      DO I = 1, K - 1
         DO II = 1, NATOM(SPC(K))
            DO JJ = 1, NATOM(SPC(I))
               RXIIJJ = RC(K,II,1) - RC(I,JJ,1)
               RYIIJJ = RC(K,II,2) - RC(I,JJ,2)
               RZIIJJ = RC(K,II,3) - RC(I,JJ,3)
 
               RIIJJSQ = RXIIJJ**2 + RYIIJJ**2 + RZIIJJ**2

               RIIJJ = DSQRT(RIIJJSQ)

               IF(RIIJJ.LE.D) GOTO 44
            END DO
         END DO
      END DO
      
      DO I = K + 1, IMLW
         DO II = 1, NATOM(SPC(K))
            DO JJ = 1, NATOM(SPC(I))
               RXIIJJ = RC(K,II,1) - RC(I,JJ,1)
               RYIIJJ = RC(K,II,2) - RC(I,JJ,2)
               RZIIJJ = RC(K,II,3) - RC(I,JJ,3)

               RIIJJSQ = RXIIJJ**2 + RYIIJJ**2 + RZIIJJ**2

               RIIJJ = DSQRT(RIIJJSQ)

               IF(RIIJJ.LE.D) GOTO 44
            END DO
         END DO
      END DO
      
      GOTO 55
      
   44 CONTINUE
   
      ACCP = .FALSE.

      WRITE(12,*)'structure bad!'
      
      RETURN
      
   55 CONTINUE
   
      ACCP = .TRUE.   

      WRITE(12,*)'structure ok!'

      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_RANGER(NOB,Rest)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      IMPLICIT none
      integer nob           ! input
      double precision rest ! output

      double precision rb, testme
      double precision NDA,NDSIG

      double precision pv, s, si
      double precision pi, smax, simax
      integer npts, i, iupper, jstar

      DIMENSION Pv(100),S(100),SI(100)

      COMMON/DIST1/RB(100),TESTME(100,2000)
      COMMON/DIST2/NDA(100),NDSIG(100)


      PI = 3.141592653589793D+00

      Rest = 0.D+00
            
C      NPTS = NOB - 1
       NPTS = NOB   
c      write(12,*) 'in dntmc_ranger npts is ',npts

c    storing distribution based on volume

      DO I = 1, NPTS
         Pv(I) = NDA(I)/(4.D+00*PI*RB(I)**2.D+00)
      ENDDO


C    smoothing distribution

      DO I = 1, 2
         CALL dntmc_SMOOTH(Pv,NPTS)
      ENDDO

c    writing distribution
c      write(12,*) 'rb and pv'
c      DO I = 1, NPTS
c         WRITE(12,*)i, RB(I),Pv(I)
c      ENDDO

c    calculating slopes via centered difference

      DO I = 2, NPTS - 1
         S(I) = (Pv(I+1) - Pv(I-1))/(RB(I+1) - RB(I-1))
      ENDDO

c     smoothing slopes
c Shawn I think this needs to be changed since s(1) is
c never set or otherwise used
c     
      s(1) = 0.0d+00
      DO I = 1, 15
          CALL dntmc_SMOOTH(S,NPTS - 2)
      ENDDO

c    writing slopes

c      write(12,*) 'rb, s'
c      DO i = 2, NPTS - 1
c         WRITE(12,*)i, RB(i),S(i)
c      ENDDO

c    finding peak in slopes

c      write(12,*) 'before dntmc_peak npts is ',npts
      CALL dntmc_PEAK(S,NPTS,SMAX,JSTAR)
      WRITE(6,*)'JSTAR =',JSTAR
      WRITE(6,*)'THE PEAK SLOPE IS AT =',RB(JSTAR),SMAX 
      WRITE(12,*)'JSTAR =',JSTAR
      WRITE(12,*)'THE PEAK SLOPE IS AT =',RB(JSTAR),SMAX 

c    refecting the slopes through the abscissa

      IUPPER = JSTAR

      WRITE(6,*)'IUPPER =',IUPPER
      WRITE(12,*)'IUPPER =',IUPPER

c      write(12,*) 'si initialized'
      DO I = 2, IUPPER
         SI(I) = -S(I)
c         write(12,*) i, si(i)
      ENDDO

c    shifting the refected slopes so that they are all positive

c      write(12,*) 'dabs of si iupper is ',iupper,dabs(si(iupper))
c      write(12,*) 'i and rb and si'
      DO I = 2, IUPPER
         SI(I) = SI(I) + DABS(SI(IUPPER))
c         WRITE(12,*)i,RB(I),SI(I)
      ENDDO

c     finding peak in positive reflected slopes

      CALL dntmc_PEAK(SI,IUPPER,SIMAX,JSTAR)
      WRITE(6,*)'JSTAR =',JSTAR
      WRITE(6,*)'THE PEAK IS AT =',RB(JSTAR),SIMAX  
      WRITE(6,*)'RSTAR =',RB(JSTAR)
      WRITE(12,*)'JSTAR =',JSTAR
      WRITE(12,*)'SIMAX =',SIMAX
      WRITE(12,*)'RSTAR =',RB(JSTAR)
        
      Rest = RB(JSTAR)

      RETURN 
      END

C***********************************************************************
      SUBROUTINE dntmc_SMOOTH(Y,NPTS)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
      double precision y  ! input/output
      integer npts        ! input
      integer imax, i
      double precision y1, ynew
      DIMENSION Y(100)

      IMAX = NPTS - 1

      Y1 = Y(1)

      DO I = 1, IMAX
         YNEW = (Y1 + 2.*Y(I) + Y(I+1))/4.
         Y1 = Y(I)
         Y(I) = YNEW
      ENDDO

      Y(NPTS) = (Y1 + 3.*Y(NPTS))/4.

      RETURN
      END

C***********************************************************************
      SUBROUTINE dntmc_PEAK(Y,NPTS,YMAX,JSTAR)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
      double precision y                ! input
      integer npts                      ! input
      double precision ymax             ! output
      integer jstar                     ! output
      integer jmax, j
      DIMENSION Y(100)

      YMAX = 0.D+00
      JMAX = NPTS - 1
   
c     write(12,*) 'npts in dntmc_peak is', npts
c     write(12,*) 'JMAX in dntmc_peak is', jmax

      DO J = 3, JMAX
         IF((Y(J).GT.Y(J-1)).AND.(Y(J).GT.YMAX))THEN
c           write(12,*) 'INSIDE J IF with Ys of ',j,y(j),y(j-1),ymax
            YMAX = Y(J)
            JSTAR = J
c           write(12,*) 'INSIDE J IF with JSTAR of ',jstar
         ENDIF
      END DO

      RETURN
      END

C***********************************************************************
C  Added NDA_OLD array
      SUBROUTINE dntmc_AVERAGEDIST(NOB,NNODES,PERCENT,nda_old)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
      double precision NDSUM,NDA,NDSIGSUM,NDSIG,NDSGMX
      double precision rb, testme
      integer nob, nnodes ! input parameters
      double precision percent     ! output
      integer i, j, jstar
C Added variables
      double precision nda_old
      dimension nda_old(100)
C End Added variables
      
      COMMON/DIST1/RB(100),TESTME(100,2000)
      COMMON/DIST2/NDA(100),NDSIG(100)


C    initialize all values to zero

c  prefer to do NOB  
c      DO J = 1, NOB - 1
       DO J = 1, NOB
         NDA(J) = 0.D+00
         NDSIG(J) = 0.D+00
C  Added initialization of nda_old
         nda_old(J) = 0.D+00
C  End added initialization
      ENDDO

C    find average for each bin across all nodes

C      DO J = 1, NOB - 1
       DO J = 1, NOB
         NDSUM = 0.D+00
         DO I = 1, NNODES
            NDSUM = NDSUM + TESTME(J,I)
         ENDDO
         NDA(J) = NDSUM/DBLE(NNODES)
      ENDDO

C    find variance for each bin across all nodes

C      DO J = 1, NOB - 1
       DO J = 1, NOB
         NDSIGSUM = 0.D+00
C  Change the variance calculation to use nda_old
C         DO I = 1, NNODES
C            NDSIGSUM = NDSIGSUM + (TESTME(J,I) - NDA(J))**2
C         ENDDO
         NDSIGSUM = (NDA(J) - NDA_OLD(J))**2
C        NDSIG(J) = DSQRT(NDSIGSUM/DBLE(NNODES))
         NDSIG(J) = NDSIGSUM
c        NDSIG(J) = DSQRT(NDSIGSUM/DBLE(NNODES - 1))
C         WRITE(6,*)RB(J),NDA(J),NDSIG(J)
      ENDDO
C  Change from maximum variance to average
         NDSIGSUM = 0.D+00
      Do J = 1, NOB
         NDSIGSUM = NDSIGSUM + NDSIG(J)
      enddo
         NDSIGSUM = DSQRT(NDSIGSUM/DBLE(NOB))
C  End changes to average variance

C  Don't use maximum error
C     FIND MAXIMUM ERROR ALONG AVERAGE DISTRIBUTION
C      
C      NDSGMX = 0.D+00
C      JSTAR = 0
C
C      DO J = 2, NOB - 1
C         IF((NDSIG(J).GT.NDSIG(J-1)).AND.(NDSIG(J).GT.NDSGMX))THEN
C            NDSGMX = NDSIG(J)
C            JSTAR = J
C         ENDIF
C      END DO
C
C      WRITE(6,*)'Maximum Error =',NDSGMX,' at J =',JSTAR
C
C      PERCENT = (NDSGMX/NDA(JSTAR))*100.D+00
c  Change Percent to an average error as a percentage
      PERCENT = NDSIGSUM*100.D+00
C  End Change

c     WRITE(6,*)'R =',RB(JSTAR),'NDA =',NDA(JSTAR),' Jstar =',JSTAR
c     WRITE(6,*)'Percent of max error at value =',PERCENT,' %'

      RETURN
      END
