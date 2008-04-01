C$Id: gibbs.F,v 1.22 2008-04-01 16:43:48 windus Exp $
C23456789012345678901234567890123456789012345678901234567890123456789012
C     modified on oct 3, 1999 for restart on mpp1
C     modified on oct 25, 1999 for random number generator restart
C     modified on may 10, 2004 for colony2a
C     modified on april 15, 2005 for ran3
C     modified on and about april 2, 2007 to clean up code (lcrosby)
C     version 5.0  Dec. 20, 2007 full release version (lcrosby)

      logical function dntmc_gibbs(rtdb)
      IMPLICIT none

C Include Statements
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "msgids.fh"
#include "msgtypesf.h"
#include "gibbsP.fh"
#include "util.fh"
#include "inp.fh"
#include "util_sgroup.fh"

C Variable declarations
      integer rtdb
      integer myid, idbig, inodesbig
      integer idmedium, idzero
      integer groupnums, g_pr
      integer iacceptt, iacceptr, iacceptv, iacceptpr, iacceptpi
      integer ivtry, ittry, irtry
      integer cont1, cont2, nmols, mprnt
      double precision temp, uu, uavr, usqrd, rsimsum, rsimavr
      logical rconf
      integer lmin, lmax
      character*4 prcfil
      character*4 prcfil2
      integer procrestart, nspc
      character*256 fprefix
      double precision rest
      logical flag
      double precision nda_old
      double precision conv_perc
      double precision rsh
      double precision loop_start, loop_time, loop_avg
      integer i,j,lnow
      double precision rminnew, rmaxnew
      external util_file_copy
      character*256 fprefixcatf
      character*256 fprefixcatt
      double precision dipole(3)
      character*256 dntmc_dir, fprefix2

C Variable declarations for common block
      integer L
      double precision rb,nda,ndsig,testme

C Common Blocks
      COMMON/MSTP/L
      COMMON/DIST1/RB(NOB),TESTME(NOB,maxcpus)
      COMMON/DIST2/NDA(NOB),NDSIG(NOB)

C Dimensions
      dimension groupnums(maxcpus)
      dimension nda_old(NOB)

C Main Program
      dntmc_gibbs = .false.
      call util_print_push
      call util_print_rtdb_load(rtdb,'dntmc')

      if(ga_nodeid() .eq. 0) then
        WRITE(6,*)'********** Starting DNTMC Module **********'
        WRITE(6,*)'*******************************************'
        call flush(6)
      endif

C Restart Post Processing
      if (rtdb_get(rtdb, 'dntmc:procrestart', MT_INT, 1,
     &    procrestart)) then
      if(ga_nodeid() .eq. 0) then
        WRITE(6,*)'********* Restart Post Processing *********'
        WRITE(6,*)'*******************************************'
        call flush(6)
      endif
        call dntmc_restart_post(procrestart,rtdb)
      goto 9898 ! Speed to end of code
      endif

C Get Info on processes
      myid = ga_nodeid()
      idbig = ga_pgroup_get_default()
      inodesbig = ga_nnodes()

C Setup GA
      if(ga_nodeid() .eq. 0) then
        WRITE(6,*)'******** Splitting Into SubGroups *********'
        WRITE(6,*)'*******************************************'
        call flush(6)
      endif

      call dntmc_setup_ga
     &(rtdb, myid, idbig, inodesbig,
     &idmedium, idzero, g_pr, groupnums) ! Now in Subgroups

C Initialize group related information
      if(ga_nodeid() .eq. 0) then
      call dntmc_setup_group
     &(rtdb, idmedium, myid, groupnums, prcfil, prcfil2, fprefix)
      else
      call dntmc_write_prcfil(idmedium, groupnums, prcfil)
      call dntmc_build_prcfil(prcfil2, myid)
      write(fprefix2,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix2))
     &     call errquit('dntmc_gibbs:rtdb get file_prefix failed'
     &     ,0,RTDB_ERR)
      write(dntmc_dir,'(256(a))') (' ', i=1,256)
      if (.not. rtdb_cget(rtdb, 'dntmc:directory',1,dntmc_dir))
     &   dntmc_dir(1:2)='./'
      write(fprefix,'(256(a))') (' ', i=1, 256)
      write(fprefix,'(3(a))')
     &     dntmc_dir(1:inp_strlen(dntmc_dir)),
     &     '/',
     &     fprefix2(1:inp_strlen(fprefix2))
      endif

C Setup Job Information
      if(ga_pgroup_nodeid(idbig) .eq. 0) then
        WRITE(6,*)'********** Setting Up DNTMC Info **********'
        WRITE(6,*)'*******************************************'
        call flush(6)
      endif

      call dntmc_setup_dntmc(rtdb, iacceptt, iacceptr,
     &  iacceptv, ivtry, ittry, irtry,
     &  iacceptpr, iacceptpi, cont1, cont2, nmols,
     &  mprnt, temp, rconf, lmin, lmax, myid, idmedium, groupnums,
     &  prcfil, uu, uavr, usqrd, rsimsum, rsimavr, nspc, dipole,
     &  conv_perc)

C Starting Monte Carlo Loop
      if(ga_pgroup_nodeid(idbig) .eq. 0) then
        WRITE(6,*)'********** Starting MC Loop ***************'
        WRITE(6,*)'*******************************************'
        call flush(6)
      endif
C      conv_perc = 0.D+00 ! Never converge 
      loop_avg = 0.00D+00
      DO i = 1, NOB
        nda_old(i) = 0.00D+00
      ENDDO
      if(util_print('Information',print_none)) then
        if(ga_pgroup_nodeid(idzero) .eq. 0) then
          WRITE(6,"(A30, ES12.3E3)")'Max MC configurations are:'
     &         ,(DBLE(lmax*util_sgroup_numgroups()))
          WRITE(6,"(A30, F12.3)")'Convergence criteria is below:'
     &         ,conv_perc
          WRITE(6,"(A7,A16, A12, A16, A16, A16, A16, A16)")
     &    '      @','      # of Total','        Rcut','       Evp. Rate'
     &    ,'      Err. Evap.','   Accept. Ratio','   Accept. Ratio',
     &     '     Convergence'
          WRITE(6,"(A7,A16, A12, A16, A16, A16, A16, A16)")
     &    '      @','  Configurations','     (Angs.)','          (s^-1)'
     &    ,'     Rate (s^-1)','       Trans (%)','         Rot (%)',
     &     '             (%)'
          WRITE(6,"(A115)")
     &'@****************************************************************
     &******************************************************'
          call flush(6)
        endif
      endif

      Do L = LMIN, LMAX
      loop_start = util_wallsec()
      call dntmc_mcloop(rtdb,temp,UU,nmols,nspc,iacceptt,iacceptr
     &  ,iacceptv,iacceptpr,iacceptpi,uavr,rsimavr,rconf,prcfil,fprefix
     &  ,mprnt,usqrd,rsimsum,idmedium,idzero,cont1,g_pr,rest,flag
     &  ,nda_old,conv_perc,ivtry,ittry,irtry,dipole)
      loop_time = util_wallsec() - loop_start
      if(loop_avg .eq. 0.00D+00) loop_avg = loop_time
      loop_avg = (loop_avg + loop_time)/2.00D+00
      if(.not. util_test_time_remaining(rtdb,NINT(loop_avg*4.00D+00)))
     & then
       if(ga_pgroup_nodeid(idbig) .eq. 0)
     &    write(6,*)'********** Not enough time remains **********'
       if(ga_pgroup_nodeid(idbig) .eq. 0) call flush(6)
       lnow = l + 1
       goto 299 ! Not Converged
      endif
      if(flag) goto 199 ! Converged
      enddo
      lnow = l
      goto 299 ! Not Converged

C End of Monte Carlo Loop 1st.
 199  CONTINUE ! Converged

C Finish up MC Loop 
      if(ga_pgroup_nodeid(idbig) .eq. 0) then
      if (util_print('fdist',print_low)) then
      WRITE(3,*)((L-1)*util_sgroup_numgroups()), 'Total Configurations'
      DO I=1, NOB
      WRITE(3,*)RB(I),NDA(I),NDSIG(I)
      ENDDO
      WRITE(3,*)'***************************************'
      call flush(3)
      endif
      endif
     
      write(fprefixcatf,'(256(a))') (' ', i=1, 256)
      write(fprefixcatf,'(3(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.restart.',
     &      prcfil(1:4)

      write(fprefixcatt,'(256(a))') (' ', i=1, 256)
      write(fprefixcatt,'(3(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.restart-1.',
     &      prcfil(1:4)

      if(ga_nodeid() .eq. 0) then
      call util_file_copy(fprefixcatf(1:inp_strlen(fprefixcatf)),
     &     fprefixcatt(1:inp_strlen(fprefixcatt)))
      endif

C Analyze Data

C Print Results

C Job Exiting Area
      goto 999 ! End Job
 299  CONTINUE ! Not Converged
C Finish up Job (Not Converged)
      IF(CONT1.EQ.1)THEN
      if (ga_nodeid() .eq. 0) then
        call dntmc_writefiles(uu, uavr, usqrd, ivtry, ittry, irtry,
     +    rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +    iacceptpi, rsh, prcfil, nmols, lnow, fprefix, dipole)
        if(util_print('debug',print_debug)) then
        write(6,*)' Restart files written for group ',prcfil(1:4)
        call flush(6)
        endif
        if(util_print('MCout',print_debug)) then
        write(12,*)' Restart files written for group ',prcfil(1:4)
        call flush(12)
        endif
      endif
      ENDIF

      if(ga_pgroup_nodeid(idbig) .eq. 0) then
      if(util_print('fdist',print_low)) then
      WRITE(3,*)(L*util_sgroup_numgroups()), 'Total Configurations'
      DO I=1, NOB
      WRITE(3,*)RB(I),NDA(I),NDSIG(I)
      ENDDO
      WRITE(3,*)'***************************************'
      call flush(3)
      endif
      endif

 999  CONTINUE ! Ending Job Correctly
C Close File Units
      if (ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_low)) then
      CLOSE(12)
      endif
      if(util_print('MCdata',print_low)) then
      CLOSE(11)
      endif
      if(util_print('Alldata',print_low)) then
      CLOSE(10)
      endif
      if (myid .eq. 0) then
      if (util_print('fdist',print_low)) then
      CLOSE(3)
      endif
      endif
      endif
C Clean up GA
      if(ga_pgroup_nodeid(idbig) .eq. 0) then
        WRITE(6,*)'********** Destroying SubGroups ***********'
        WRITE(6,*)'*******************************************'
        call flush(6)
      endif

      if (ga_pgroup_nodeid(idmedium) .eq. 0) then
      call ga_pgroup_set_default(idzero)
      if(.not. ga_destroy(g_pr)) then
        call errquit('Failure to clean up GA',0,GA_ERR)
      endif
      endif
C Finalize Subgroups
      call ga_pgroup_set_default(idmedium)
      call util_sgend(rtdb)
      call util_print_pop
C End Main Program
 9898 dntmc_gibbs = .true.
      return
      end
C***********************************************************************

C************************************************************************
      Subroutine dntmc_writefiles(uu, uavr, usqrd, ivtry, ittry, irtry,
     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, fprefix, dipole)
C***********************************************************************
      Implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "inp.fh"

      double precision uu, usum, uavr, usqrdsum, usqrd
      double precision rsim, rsimsum, rsimavr
      double precision seed, r, numbersum, numberavr
      double precision rmin, rmax, rinc, rcm, rs, dipole(3)
      integer number, iacceptt, iacceptv, iacceptr, iacceptpr, iacceptpi
      integer ivtry, ittry, irtry
      integer rsmcstps, i, j, ir
      integer natom, spc
      double precision rsh
      integer nmols
      character*4 prcfil
      character*256 fprefix
      character*256 fprefixcat
      character*16 tag
      integer ti

      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/AVR/USUM,USQRDSUM
      COMMON/NUMBER1/NUMBER(NOB)
      COMMON/NUMBER2/NUMBERSUM(NOB),NUMBERAVR(NOB)
      COMMON/SEEDS/SEED(1)
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)
      COMMON/ATTAG/tag(mxatom*mxmols)


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
         WRITE(20,*)(DIPOLE(IR), IR=1,3)
         WRITE(20,*)RMIN, RMAX
         WRITE(20,*)RSIM
         WRITE(20,*)RSIMSUM
         WRITE(20,*)RSIMAVR
         WRITE(20,*)IACCEPTT
         WRITE(20,*)IACCEPTR
         WRITE(20,*)IACCEPTV
         WRITE(20,*)IVTRY,ITTRY,IRTRY
         WRITE(20,*)IACCEPTPR
         WRITE(20,*)IACCEPTPI
         WRITE(20,*)rsmcstps
         ti = 1
         DO I = 1, NMOLS
            DO J = 1, NATOM(SPC(I))
               WRITE(20,*)tag(ti),(R(I,J,IR), IR = 1,3)
               ti = ti + 1
            ENDDO
         ENDDO

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
C  End of Subroutine dntmc_writefiles(uu, uavr, usqrd, ivtry, ittry, irtry,
C     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
C     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, fprefix, dipole)
C***********************************************************************

C***********************************************************************
      Subroutine dntmc_mcloop(rtdb,temp,UU,nmols,nspc,iacceptt,iacceptr
     &  ,iacceptv,iacceptpr,iacceptpi,uavr,rsimavr,rconf,prcfil,fprefix
     &  ,mprnt,usqrd,rsimsum,idmedium,idzero,cont1,g_pr,rest,flag
     &  ,nda_old,conv_perc,ivtry,ittry,irtry,dipole)
C***********************************************************************
      Implicit none
C Include Statements
#include "gibbsP.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "msgtypesf.h"
#include "util.fh"
#include "rtdb.fh"
#include "util_sgroup.fh"

C Variable Declarations
      double precision temp  !input
      logical acceptt, acceptr, acceptv, acceptpr, acceptpi ! internal
      integer rtdb, nmols !input
      integer iacceptt, iacceptr, iacceptv, iacceptpr, iacceptpi !input/output
      integer ivtry, ittry, irtry ! input/output
      double precision UU, UAVR, CV, RSIMAVR !input/output
      integer nspc !input
      logical rconf ! input
      double precision rinc, rsh !internal
      integer lnow ! internal
      integer mprnt ! input
      character*4 prcfil !input
      character*256 fprefix ! input
      double precision usqrd, rsimsum ! input/output
      integer idmedium, idzero, cont1 ! input
      integer ilo,ihi,jlo,jhi,ld(2) ! internal
      integer g_pr !input
      double precision percent ! internal
      double precision rest !output
      logical flag !output
      double precision aratiott, aratiotr, aratiotv, aratiotpr,aratiotpi
      double precision aratiovtry, aratiottry, aratiortry
      integer i,j,s,t, z
      double precision conv_perc ! input
      double precision nda_old ! input/output
      dimension nda_old(NOB) ! input/output
      double precision dipole(3) !input/output
      integer movetry
      double precision sumx
      double precision flux(NOB), flux2(NOB)
      double precision ERRX
      double precision sumerr(NOB), sumerr2(NOB)
      double precision press
      logical inflex
      double precision rflux(NOB), ratec(NOB)
      double precision errrate, intpress, ratef, rbf
      double precision acceptratio
      logical dynamic

C Variables for common blocks
      integer L, ncurr,nmin,nmax
      double precision rmin, rmax, rsim
      double precision RB, TESTME
      double precision numbersum,numberavr
      double precision nda, ndsig
      double precision e,k,esqrd,kb,PI
      double precision m, mtot
      integer uspc
      double precision disp, alpha
      integer rotf

C Common Blocks
      COMMON/MSTP/L
      COMMON/GRAND/NMIN,NCURR,NMAX
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/DIST1/RB(NOB),TESTME(NOB,maxcpus)
      COMMON/NUMBER2/NUMBERSUM(NOB),NUMBERAVR(NOB)
      COMMON/DIST2/NDA(NOB),NDSIG(NOB)
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/MASSES/M(MXSP,MXATOM),MTOT(MXSP)
      COMMON/USPECIES/USPC(MXSP)
      COMMON/TRPS/DISP,ALPHA,ROTF

C Main Program
         press = 1.7D-05
         ACCEPTT = .FALSE. ! translate
         ACCEPTR = .FALSE. ! rotation
         ACCEPTV = .FALSE. ! volume
         ACCEPTPR = .FALSE. !????
         ACCEPTPI = .FALSE.
         movetry = 0
          if (.not. rtdb_get(rtdb, 'dntmc:dynamicsteps', mt_log, 1
     &       , dynamic))
     $       dynamic = .false.

C do one move, calculate energy and return the acceptance ratios and energy
      if(rconf) then
         CALL dntmc_MOVE_RCONF(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
     +             ,ACCEPTPI,UU, rtdb, nmols, movetry, dipole)
      else
         CALL dntmc_MOVE(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
     +             ,ACCEPTPI,UU, rtdb, nmols, movetry)
      endif
         CALL dntmc_ROOTMEANSQUARE(L,NCURR,NSPC)

         IF(ACCEPTT) IACCEPTT = IACCEPTT + 1
         IF(ACCEPTR) IACCEPTR = IACCEPTR + 1
         IF(ACCEPTV) IACCEPTV = IACCEPTV + 1
         IF(ACCEPTPR) IACCEPTPR = IACCEPTPR + 1
         IF(ACCEPTPI) IACCEPTPI = IACCEPTPI + 1

         IF(movetry .eq. 1) then
           IVTRY = IVTRY + 1
         ELSE IF(movetry .eq. 2) then
           ITTRY = ITTRY + 1
         ELSE IF(movetry .eq. 4) then
           IRTRY = IRTRY + 1
         ELSE IF(movetry .eq. 3) then
           IVTRY = IVTRY + 1
           ITTRY = ITTRY + 1
         ELSE IF(movetry .eq. 5) then
           IVTRY = IVTRY + 1
           IRTRY = IRTRY + 1
         ELSE IF(movetry .eq. 6) then
           ITTRY = ITTRY + 1
           IRTRY = IRTRY + 1
         ELSE IF(movetry .eq. 7) then
           IVTRY = IVTRY + 1
           ITTRY = ITTRY + 1
           IRTRY = IRTRY + 1
         ENDIF

         ARATIOTT = (DBLE(IACCEPTT)/DBLE(ITTRY))*100.D+00 
         ARATIOTR = (DBLE(IACCEPTR)/DBLE(IRTRY))*100.D+00 
         ARATIOTV = (DBLE(IACCEPTV)/DBLE(IVTRY))*100.D+00
         IF (ITTRY .eq. 0) ARATIOTT = 0.D+00
         IF (IRTRY .eq. 0) ARATIOTR = 0.D+00
         IF (IVTRY .eq. 0) ARATIOTV = 0.D+00
         ARATIOVTRY = (DBLE(IVTRY)/DBLE(L))*100.D+00
         ARATIOTTRY = (DBLE(ITTRY)/DBLE(L))*100.D+00
         ARATIORTRY = (DBLE(IRTRY)/DBLE(L))*100.D+00
          
         ARATIOTPR = (DBLE(IACCEPTPR)/DBLE(L))*100.D+00
         ARATIOTPI = (DBLE(IACCEPTPI)/DBLE(L))*100.D+00

         IF(rconf) then
         ACCEPTRATIO=(1.0D+00/3.0D+00)*
     &               (ARATIOTT + ARATIOTR + 100.0D+00)
         else
         ACCEPTRATIO=(1.0D+00/4.0D+00)*
     &               (ARATIOTT+ARATIOTR+ARATIOTV+100.0D+00)
         endif

C calculate average of stuff
         CALL dntmc_AVERAGE(L,TEMP,UU,UAVR,Cv,RSIMAVR,usqrd,rsimsum)
c calculate distribution of different radii of the sphere
         CALL dntmc_NUMBERHIST(L,rconf)

       RINC = (RMAX - RMIN)/DBLE(NOB)
       RSH = RMIN
       DO I = 1, NOB
          RSH = RSH + RINC
          RB(I) = RSH
       ENDDO

C Calculate Radial Distribution
         CALL dntmc_NUMBERAVERAGE(L,temp,rconf)
         flag = .false.
C Occasional grouping
      IF (MOD(L,MPRNT).EQ.0)THEN
       if(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_debug)) then
      WRITE(12,*)'*******************************'
      WRITE(12,*)'THIS IS THE END OF STEP = ',L
      WRITE(12,*)'the number of molecules currently'
      WRITE(12,*)'in the simulation are',NCURR
      WRITE(12,*)'ARATIOT = ',ARATIOTT,' %'
      WRITE(12,*)'ARATIOR = ',ARATIOTR,' %'
      WRITE(12,*)'ARATIOV = ',ARATIOTV,' %'
      WRITE(12,*)'ARATIOVTRY = ',ARATIOVTRY,' %'
      WRITE(12,*)'ARATIOTTRY = ',ARATIOTTRY,' %'
      WRITE(12,*)'ARATIORTRY = ',ARATIORTRY,' %'
      WRITE(12,*)'ARATIOPR = ',ARATIOTPR,' %'
      WRITE(12,*)'ARATIOPI = ',ARATIOTPI,' %'
      WRITE(12,*)'UNEW =',UU,' kcal/mole '
      WRITE(12,*)'RSIM =',RSIM,' ANGSTROMS'
      WRITE(12,*)'<RSIM> =',RSIMAVR,' ANGSTROMS'
      WRITE(12,*)'<U> =',UAVR,' kcal/mole'
      WRITE(12,*)'<Cv> =',Cv,' kcal/mole-K'
      WRITE(12,*)'*******************************'
      call flush(12)
      endif
       endif

      IF (Dynamic) then
      IF (ARATIOTT .lt. 75.0D+00) then
        DISP = DISP*0.90D+00
      Else if (ARATIOTT .gt. 75.0D+00) then
        DISP = DISP*1.10D+00
      endif

      IF (ARATIOTR .lt. 75.0D+00) then
        ALPHA = ALPHA*0.90D+00
      else if (ARATIOTR .gt. 75.0D+00) then
        ALPHA = ALPHA*1.10D+00
      endif
      endif

      IF (Mod(L,MPRNT*NOB) .eq. 0) then
C Grouping Together Group Information
      if (ga_nodeid().eq.0) then
         call ga_pgroup_set_default(idzero)
         ilo = 1
         ihi = NOB
         jlo = ga_nodeid() + 1
         jhi = ga_nodeid() + 1

         ld(1)  = NOB   !must be physical dimension of local array
         ld(2)  = 1
         call ga_put(g_pr,ilo,ihi,jlo,jhi,NUMBERAVR,ld)
         call ga_sync
         jlo = 1
         jhi = ga_nnodes()
         call ga_get(g_pr,ilo,ihi,jlo,jhi,TESTME,ld)
         CALL dntmc_AVERAGEDIST(ga_nnodes(),PERCENT,nda_old)
      endif  ! idzero section
      call ga_pgroup_set_default(idmedium)
      call ga_brdcst(msg_dntmc0+MSGDBL, PERCENT, mdtob(1), 0)
      call ga_brdcst(msg_dntmc1+MSGDBL, NDA, mdtob(NOB), 0)
      call ga_brdcst(msg_dntmc2+MSGDBL, NDSIG, mdtob(NOB), 0)

      if(ga_pgroup_nodeid(idzero) .eq. 0) then
      if(util_print('fdist',print_low)) then
      WRITE(3,*)(L*util_sgroup_numgroups()), 'Total Configurations'
      endif
      if(util_print('fdist',print_low)) then
      DO I=1, NOB
      WRITE(3,*)RB(I),NDA(I),NDSIG(I)
      ENDDO
      WRITE(3,*)'***************************************'
      call flush(3)
      endif
C      if(util_print('information',print_none)) then
C      WRITE(6,*)''
C      WRITE(6,*)(L*util_sgroup_numgroups()), 'Total Configurations Done'
C      call flush(6)
C      endif
      endif

C  Added save of old NDA array
         Do i=1, NOB
         nda_old(i) = NDA(i)
         enddo
C  End added old save

C      if(ga_pgroup_nodeid(idzero) .eq. 0) then
C         if(util_print('Information',print_none)) then
C         write(6,*)'Max percent difference in Markov Chains is',percent
C         call flush(6)
C         endif
C      endif
      IF(PERCENT .LT. conv_perc) then
         flag = .true.
          if(ga_pgroup_nodeid(idzero) .eq. 0) then
          if(util_print('Information',print_none)) then
          write(6,'(a)')'Distribution has converged below', conv_perc
          call flush(6)
          endif
          endif
      Else
         flag = .false.
C         if(ga_pgroup_nodeid(idzero) .eq. 0) then
C         if(util_print('Information',print_none)) then
C         write(6,*)'Distribution has not converged below', conv_perc
C         call flush(6)
C         endif
C         endif
      endif   

C      if(ga_nodeid() .eq. 0) write(6,*)
C     &   'percent is ', percent, ga_nodeid(), idmedium
C      IF(PERCENT.LT.conv_perc)THEN
C       if(ga_nodeid() .eq. 0) then
C          WRITE(6,*)'the average distribution has converged', idmedium
C          WRITE(6,*)'PERCENT =',PERCENT,idmedium
C          call flush(6)
C          if(util_print('MCout',print_debug)) then
C          WRITE(12,*)'the average distribution has converged'
C          WRITE(12,*)'PERCENT =',PERCENT
C          call flush(12)
C          endif
C       endif
Cc gather best estimate of R
C          CALL dntmc_RANGER(Rest, rconf)
C        if(ga_nodeid() .eq. 0) then
C          WRITE(6,*)'best estimate of Rcut =',Rest,idmedium
C          call flush(6)
C          if(util_print('MCout',print_debug)) then
C          WRITE(12,*)'best estimate of Rcut =',Rest
C          call flush(12)
C          endif
C        endif
C        flag = .true.
C        
C      ELSE
C       if(ga_nodeid() .eq. 0) then
C          WRITE(6,*)'the average distribution has not converged yet'
C     &    , idmedium
C          WRITE(6,*)'PERCENT =',PERCENT,idmedium
C          call flush(6)
C          if(util_print('MCout',print_debug)) then
C          WRITE(12,*)'the average distribution has not converged yet'
C          WRITE(12,*)'PERCENT =',PERCENT
C          call flush(12)
C          endif
C       endif
C
Cc gather best estimate of R
C          CALL dntmc_RANGER(Rest, rconf)
C        if(ga_nodeid() .eq. 0) then
C          WRITE(6,*)'best estimate of Rcut =',Rest,idmedium
C          call flush(6)
C        if(util_print('MCout',print_debug)) then
C          WRITE(12,*)'best estimate of Rcut =',Rest
C          call flush(12)
C        endif
C        endif
C        flag = .false.
C
C      ENDIF

C  Rconfig Data Analysis
      IF (rconf) then
C  Calculate Flux
        flux(1) = NDA(1)
        Do i=2, NOB
          flux(i) = NDA(i) + flux(i-1)
        End Do ! CPF in flux
        Do i = 1, NOB
          flux2(i) = NDA(i)
        End Do ! Derivative CPF in flux2
C  Errors in Flux
        sumerr(1) = NDSIG(1) * NDSIG(1)
        Do i = 2, NOB
          sumerr(i) = (NDSIG(i)*NDSIG(i)) + sumerr(i-1)
        END DO
        Do i = 1, NOB
          sumerr(i) = DSQRT(sumerr(i))
        End DO !Error in CPF in sumerr
        Do i = 1, NOB
          sumerr2(i) = NDSIG(i)
        ENDDO !Error in Derivative CPF in sumerr2
      ELSE
        Do i = 1, NOB
          flux(i) = NDA(i)/(4.D+00*PI*(RB(i)**2.D+00))
          flux(i) = flux(i)*DEXP((4.D+00/3.D+00*PI*(RB(i)**3.D+00))
     &              *press/kb/temp)
        ENDDO ! CPF in flux
        Do i = 1, NOB
          if ((i .eq. 1) .or. (i .eq. NOB)) then
            flux2(i) = 0.00D+00
          else
          flux2(i) = (flux(i+1) - flux(i-1))/2.00D+00
          endif
        ENDDO ! Derivative CPF in flux2
C  Errors in Flux
        Do i = 1, NOB
          sumerr(i) = NDSIG(i)/(4.D+00*PI*(RB(i)**2.D+00))
          sumerr(i) = sumerr(i)*DEXP((4.D+00/3.D+00*PI*(RB(i)**3.D+00))
     &              *press/kb/temp)
        ENDDO !Error in CPF in sumerr
        DO i = 1, NOB
          if ((i .eq. 1) .or. (i .eq. NOB)) then
            sumerr2(i) = 0.00D+00
          else
            sumerr2(i) =(sumerr(i-1)**2.00D+00)+(sumerr(i+1)**2.00D+00)
            sumerr2(i) = DSQRT(sumerr2(i))/2.00D+00
          endif
        ENDDO !Error in Derivative CPF in sumerr2
      ENDIF
C  Calculate ERROR
      ERRX = 0.00D+00
      Do i=1, NOB - 1
        ERRX = ERRX + RB(i+1) - RB(i)
      End Do
      ERRX = ERRX/DBLE(NOB - 1) ! Error related to bin interval
C  Reactive Flux
      Do i=1, NOB
        rflux(i) = flux2(i)*kb*temp/ERRX
      enddo ! Reactive flux in rflux
C Rate Constants
      Do i=1, NOB
        if(flux(i) .eq. 0.D+00) then
        ratec(i) = 1.D+99
        else
        ratec(i) = 1.D+00*rflux(i)/flux(i)*1.D+13/
     &  DSQRT(2.D+00*PI*kb*temp*mtot(uspc(1))/4.184D+00) ! only first species
        endif
      Enddo ! Rate Constants in ratec
C Find Minimum Flux        From Front to Back
C Minima is point i .lt. i-1 .and. .le. i+1 (i+1 .gt. 0)
C Global minimum is point i .lt. sumx
C Find Maximum Flux
      sumx = rflux(1)
      j = 1
      Do i = 2, NOB
        IF (rflux(i) .gt. sumx) then
          sumx = rflux(i)
          j = i
        EndIf
      EndDo
      if (j .eq. NOB) j = NOB -1
      if (j .eq. 1) j = 2
      z = j

C  Find Minimum Flux
      j = z
      if (j .gt. NOB/2) then

      DO i = z , 2, -1
        IF ( (rflux(i) .lt. rflux(i-1)) .and.
     &       (rflux(i) .le. rflux(i+1)) .and.
     &       (rflux(i) .lt. sumx)       .and.
     &       (rflux(i+1) .gt. 0.00D+00) ) then
          sumx = rflux(i)
          j = i !indux in j
        endif
      End Do

      else

      DO i = z , NOB - 1
        IF ( (rflux(i) .lt. rflux(i-1)) .and.
     &       (rflux(i) .le. rflux(i+1)) .and.
     &       (rflux(i) .lt. sumx)       .and.
     &       (rflux(i+1) .gt. 0.00D+00) ) then
          sumx = rflux(i)
          j = i ! index in j
        endif
      End Do

      endif
C Results
      errrate = ratec(j)*DSQRT((sumerr(j)/flux(j))**2.D+00
     &          + (flux2(j)/sumerr2(j))**2.D+00)
      intpress = rflux(j)/flux(j)/4/PI/(RB(j)**2.00D+00)
      intpress = intpress*kb*temp
      intpress = intpress*4.184*1.00D+30/1.01325D+02/6.022142D+23
      ratef = ratec(j)
      rbf = RB(j)

C Print Results
      if(util_print('Information',print_none)) then
        if(ga_pgroup_nodeid(idzero) .eq. 0) then
          WRITE(6,"(A7,I16, F12.3, ES16.3E3, ES16.3E3, F16.3, F16.3,
     &    F16.3)")
     &    '   @   ',(L*util_sgroup_numgroups()), rbf, ratef, errrate
     &    ,aratiott,aratiotr, percent
          call flush(6)
        endif
      endif
C      if(ga_nodeid() .eq. 0) then
C      WRITE(6,"(A, I5, A, ES18.10E3, A, ES18.10E3, I5)")
C     &  'best estimate of evap. rate (flux) for species ',i,' = '
C     &  ,sumx,' with error = ',sumy,idmedium
C      call flush(6)

CC Calculate Flux
C      IF (rconf) then
C        flux(1) = NDA(1)
C        Do i=2, NOB
C          flux(i) = NDA(i) + flux(i-1)
C        End Do
C        Do i=1, NOB
C          if (flux(i) .eq. 0.D+00) then
C            flux(i) = 1.D+10
C          else
C            flux(i) = -1.D+00*kb*temp*DLOG(flux(i))
C          endif
C        end do ! Now a Helmholtz Free Energy in flux
C      ELSE
C        Do i=1, NOB
C         flux(i) = NDA(i)/(4.D+00*PI*(RB(i)**2.D+00))
C         if (flux(i) .eq. 0.D+00) then
C           flux(i) = 1.D+10
C         else
C           flux(i) = -1.D+00*kb*temp*DLOG(flux(i))
C         endif
C         flux(i) = flux(i) - (4.D+00/3.D+00*PI*(RB(i)**3.D+00)*press)
C        END Do   ! Now a Helmloltz Free Energy in flux
C      ENDIF
CC Calculate Error
C      ERRX = 0.00D+00
C      Do i=1, NOB - 1
C        ERRX = ERRX + RB(i+1) - RB(i)
C      End Do
C      ERRX = ERRX/DBLE(NOB - 1) ! Error related to bin interval
C      ERRX = DABS(ERRX)
C 
C      IF (rconf) then
C        sumnda(1) = NDA(1)
C        sumerr(1) = NDSIG(1)*NDSIG(1)
C        Do i = 2, NOB
C          sumnda(i) = NDA(i) + sumnda(i-1)
C          sumerr(i) = (NDSIG(i)*NDSIG(i)) + sumerr(i-1)
C        END DO
C        Do i = 1, NOB
C          sumerr(i) = DSQRT(sumerr(i))
C        End DO
C        Do i = 1, NOB
C          if (sumnda(i) .eq. 0.D+00) then
C            ERRF1(i) = 0.D+00
C          else
C            ERRF1(i) = DABS(kb*temp*sumerr(i)/sumnda(i))
C          endif
C        END DO ! have error in Helmholtz Free Energy in ERRF1
C      ELSE
C        Do i = 1, NOB
C          if (NDA(i) .eq. 0.D+00) then
C            ERRF1(i) = 0.D+00
C          else
C            ERRF1(i) = ((NDSIG(i)/NDA(i))**2.D+00)
C          endif
C          if(RB(i) .eq. 0.D+00) then
C            ERRF1(i) = ERRF1(i) + 0.D+00
C          else
C            ERRF1(i) = ERRF1(i) + (4.D+00*((ERRX/RB(i))**2.D+00))
C          endif
C          ERRF1(i) = DSQRT(ERRF1(i))
C          ERRF1(i) = ERRF1(i)*kb*temp
C          ERRF1(i) = ERRF1(i)**2.D+00
C          if (RB(i) .eq. 0.D+00) then
C            ERRF1(i) = ERRF1(i) + 0.D+00
C          else
C            ERRF1(i) = ERRF1(i) + (9.D+00*((press*4.D+00/3.D+00*PI
C     &      *(RB(i)**3.D+00))**2.D+00)
C     &      *((ERRX/RB(i))**2.D+00))
C          endif
C          ERRF1(i) = DSQRT(ERRF1(i))
C        End Do ! have error in Helmholtz Free Energy in ERRF1
C      Endif
C
CC Loop to reduce grid points
C      Do s=1, (NOB/10)
C         j = 1
C         Do i=1, NOB, s
C           flux2(j) = 0.D+00
C           RB2(j) = 0.D+00
C           ERRF2(j) = 0.D+00
C           Do t=i, i+s-1
C             flux2(j) = flux2(j) + flux(t)
C             RB2(j) = RB2(j) + RB(t)
C             ERRF2(j) = ERRF2(j) + (ERRF1(t)*ERRF1(t))
C           End do
C           flux2(j) = flux2(j)/DBLE(s)
C           RB2(j) = RB2(j)/DBLE(s)
C           ERRF2(j) = DSQRT(ERRF2(j))/DBLE(s)
C           ERRX2 = ERRX/DSQRT(DBLE(s))
C           j = j + 1
C         End do ! information reduced
CC New bin definitions in RB2, Helholtz Free Energy in flux2, and error in ERRF2
CC New bin error in ERRX2
C            
C      Do i=1,NOB
C        slope(i) = 0.D+00
C        ERRS(i) = 0.D+00
C      End do
C      Do i = 2, (NOB/s) - 1
C       sumx = 0.D+00
C       sumy = 0.D+00
C       sumxy = 0.D+00
C       sumx2 = 0.D+00
C       sumydx = 0.D+00
C       sumxdy = 0.D+00
C       sumdx = 0.D+00
C       sumdy = 0.D+00
C       sumxdx = 0.D+00
C       Do j = i-1, i+1
C         sumx = sumx + RB2(j)
C         sumx2 = sumx2 + (RB2(j)*RB2(j))
C         sumy = sumy + flux2(j)
C         sumxy = sumxy + (RB2(j)*flux2(j))
C       End Do
C      slope(i) = ((3.D+00*sumxy) - (sumx*sumy))/
C     &           ((3.D+00*sumx2) - (sumx*sumx))   ! Have derivative of A in slope
C       Do j = i-1, i+1
C         sumdx = sumdx + ((ERRX2*((3.D+00*flux2(j))-sumy-
C     &          (slope(i)*((6.D+00*RB2(j))-(2.D+00*sumx)))))**2.D+00)
C         sumdy = sumdy + ((ERRF2(j)*((3.D+00*RB2(j))-sumx))**2.D+00)
C       End Do
C      ERRS(i) = DSQRT(sumdx + sumdy)/((3.D+00*sumxy)-(sumx*sumy))
C      ERRS(i) = DABS(slope(i)*ERRS(i))
C      End Do
C
C      if(rconf) then
C      Do i=1,NOB
C        slope2(i) = 0.D+00
C      End do
C      Do i = 3, (NOB/s) - 2
C       sumx = 0.D+00
C       sumy = 0.D+00
C       sumxy = 0.D+00
C       sumx2 = 0.D+00
C       Do j = i-1, i+1
C         sumx = sumx + RB2(j)
C         sumx2 = sumx2 + (RB2(j)*RB2(j))
C         sumy = sumy + slope(j)
C         sumxy = sumxy + (RB2(j)*slope(j))
C       End Do
C      slope2(i) = ((3.D+00*sumxy) - (sumx*sumy))/
C     &           ((3.D+00*sumx2) - (sumx*sumx))   ! Have second derivative of A in slope2
C      End Do
C
C      inflex = .false.
C      Do i=3, (NOB/s) - 2
C        if ((slope2(i) .gt. 0.D+00) .and.
C     &     (inflex .eq. .true.)) goto 987 ! do not exit loop
C        if (slope2(i) .lt. 0.D+00) inflex = .true.
C      End Do
C      Else
C      Do i=2, (NOB/s) - 1
C        if (slope(i) .gt. 0.D+00) goto 987 ! do not exit loop
C      End Do
C      ENDIF
C      GOTO 999 ! Exit Loop
C 987  CONTINUE
C      End Do ! Reduce Grid Loop
C      s = s - 1
C 999  CONTINUE
C   
C      DO i = 2, (NOB/s) - 1
C        flux2(i) = -1.D+00*DEXP(-1.D+00*flux2(i)/(kb*temp))*slope(i) ! Now a Reactive Flux
C      END DO
C     
C      if(ga_nodeid() .eq. 0) then
C      WRITE(6,*)'Grid points have been reduced from', NOB, 'to',NOB/s
C     & , idmedium
C      call flush(6)
C      if(util_print('MCout',print_debug)) then
C      WRITE(12,*)'Grid points have been reduced from', NOB, 'to',NOB/s
C      call flush(12)
C      endif
C      endif
C
CC Find Minimum Flux
C      sumx = flux2(2) ! using sumx to hold min flux
C      j = 2
C      DO i = 3 , (NOB/s) -2
C        IF (((flux2(i) .lt. flux2(i-1)).and.
C     &    ((flux2(i) .lt. sumx)).or.(sumx .eq. 0.00D+00))
C     &    .and. (flux2(i) .le. flux2(i+1))) then
C          sumx = flux2(i)
C          j = i
C        endif
C      End Do
C      if(ga_nodeid() .eq. 0) then
C      WRITE(6,"(A,F16.6,I5)")
C     & 'best estimate of Rcut (flux)= ',RB2(j),idmedium
C      call flush(6)
C      if(util_print('MCout',print_debug)) then
C      WRITE(12,"(A,F16.6)")'best estimate of Rcut (flux)= ',RB2(j)
C      call flush(12)
C      endif
C      endif
C
CC evaporation rate constant rate constant
CC Using sumx to hold rate
CC Using sumy to hold error
C      Do i=1, nspc
C      sumx = -1.D+13*slope(j)/DSQRT(2.D+00*PI*kb*temp*mtot(uspc(i))
C     & /4.184D+00)
C      sumy = sumx/slope(j)*ERRS(j)
C      sumy = DABS(sumy)
C
C      if(ga_nodeid() .eq. 0) then
C      WRITE(6,"(A, I5, A, ES18.10E3, A, ES18.10E3, I5)")
C     &  'best estimate of evap. rate (flux) for species ',i,' = '
C     &  ,sumx,' with error = ',sumy,idmedium
C      call flush(6)
C      if(util_print('MCout',print_debug)) then
C      WRITE(12,"(A, I5, A, ES18.10E3, A, ES18.10E3)")
C     &  'best estimate of evap. rate (flux) for species ',i,' = '
C     &  ,sumx,' with error = ',sumy
C      call flush(12)
C      endif
C      endif
C      End Do
      ENDIF ! ON MOD(L,mprnt*nob)

C Restart files
      lnow = L + 1
      IF(CONT1.EQ.1)THEN
      if (ga_nodeid() .eq. 0) then
        call dntmc_writefiles(uu, uavr, usqrd, ivtry, ittry, irtry,
     +    rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +    iacceptpi, rsh, prcfil, nmols, lnow, fprefix, dipole)
        if(util_print('debug',print_debug)) then
        write(6,*)' Restart files written for group ',prcfil(1:4)
        call flush(6)
        endif
        if(util_print('MCout',print_debug)) then
        write(12,*)' Restart files written for group ',prcfil(1:4)
        call flush(12)
        endif
      endif
      ENDIF

      ENDIF ! ON MOD(L,mprnt)
C End Subroutine
      return
      end
C***********************************************************************
C End of Subroutine dntmc_mcloop(rtdb,temp,UU,nmols,nspc,iacceptt,iacceptr
C     &  ,iacceptv,iacceptpr,iacceptpi,uavr,rsimavr,rconf,prcfil,fprefix
C     &  ,mprnt,usqrd,rsimsum,idmedium,idzero,cont1,g_pr,rest,flag
C     &  ,nda_old,conv_perc,ivtry,ittry,irtry,dipole)
C***********************************************************************
C***********************************************************************
      SUBROUTINE dntmc_RANGER(Rest, rconf)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      IMPLICIT none
#include "gibbsP.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "util_sgroup.fh"

      double precision rest ! output

      double precision rb, testme
      double precision NDA,NDSIG

      double precision pv, s, si
      double precision pi, smax, simax
      integer npts, i, iupper, jstar
      logical rconf
      double precision hold

      DIMENSION Pv(NOB),S(NOB),SI(NOB)

      COMMON/DIST1/RB(NOB),TESTME(NOB,maxcpus)
      COMMON/DIST2/NDA(NOB),NDSIG(NOB)


      PI = 3.141592653589793D+00

      Rest = 0.D+00

C      NPTS = NOB - 1
       NPTS = NOB
c      write(12,*) 'in dntmc_ranger npts is ',npts

c    storing distribution based on volume
      IF (rconf) then
        Pv(1) = NDA(1)
        DO I = 2, NOB
          Pv(I) = NDA(I) + Pv(I-1)
        End Do
      Else
        DO I = 1, NPTS
          Pv(I) = NDA(I)/(4.D+00*PI*(RB(I)**2.D+00))
        ENDDO
        hold = Pv(1)
        DO I = 2, NOB
          if (Pv(I) .ge. hold) then
            hold = Pv(I)
          endif
        End Do
        DO I = 1, NOB
          Pv(I) = Pv(I)/hold
        End Do
      Endif          


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
      s(NPTS) = 0.0D+00
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
      if(ga_nodeid() .eq. 0) then
      WRITE(6,*)'JSTAR =',JSTAR
      WRITE(6,*)'THE PEAK SLOPE IS AT =',RB(JSTAR),SMAX
      call flush(6)
      if(util_print('MCout',print_debug)) then
      WRITE(12,*)'JSTAR =',JSTAR
      WRITE(12,*)'THE PEAK SLOPE IS AT =',RB(JSTAR),SMAX
      call flush(12)
      endif
      endif

c    refecting the slopes through the abscissa

      IUPPER = JSTAR

      if(ga_nodeid() .eq. 0) WRITE(6,*)'IUPPER =',IUPPER
      if(ga_nodeid() .eq. 0) call flush(6)
      if(util_print('MCout',print_debug)) then
      if(ga_nodeid() .eq. 0) WRITE(12,*)'IUPPER =',IUPPER
      if(ga_nodeid() .eq. 0) call flush(12)
      endif


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
      if(ga_nodeid() .eq. 0) then
      WRITE(6,*)'JSTAR =',JSTAR
      WRITE(6,*)'THE PEAK IS AT =',RB(JSTAR),SIMAX
      WRITE(6,*)'RSTAR =',RB(JSTAR)
      call flush(6)
      if(util_print('MCout',print_debug)) then
      WRITE(12,*)'JSTAR =',JSTAR
      WRITE(12,*)'SIMAX =',SIMAX
      WRITE(12,*)'RSTAR =',RB(JSTAR)
      call flush(12)
      endif
      endif
      Rest = RB(JSTAR)

      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_RANGER(Rest, rconf)
C***********************************************************************

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
      DIMENSION Y(NPTS)

      YMAX = 0.D+00
      JMAX = NPTS - 1
      JSTAR = 3

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
C End of SUBROUTINE dntmc_PEAK(Y,NPTS,YMAX,JSTAR)
C***********************************************************************

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
      DIMENSION Y(NPTS)

      IMAX = NPTS - 1

      Y1 = Y(1)

      DO I = 1, IMAX
         YNEW = (Y1 + 2.D+00*Y(I) + Y(I+1))/4.D+00
         Y1 = Y(I)
         Y(I) = YNEW
      ENDDO

      Y(NPTS) = (Y1 + 3.D+00*Y(NPTS))/4.D+00

      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_SMOOTH(Y,NPTS)
C***********************************************************************

C***********************************************************************
C  Added NDA_OLD array
      SUBROUTINE dntmc_AVERAGEDIST(NNODES,PERCENT,nda_old)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"
#include "util_sgroup.fh"


      double precision NDSUM,NDA,NDSIGSUM,NDSIG,NDSGMX
      double precision rb, testme
      integer nnodes ! input parameters
      double precision percent     ! output
      integer i, j, jstar
C Added variables
      double precision nda_old
      dimension nda_old(NOB)
      double precision meandist(maxcpus+1)
      double precision sumdist(maxcpus+1)
C End Added variables

      COMMON/DIST1/RB(NOB),TESTME(NOB,maxcpus)
      COMMON/DIST2/NDA(NOB),NDSIG(NOB)


C    initialize all values to zero

c  prefer to do NOB
c      DO J = 1, NOB - 1
       DO J = 1, NOB
         NDA(J) = 0.D+00
         NDSIG(J) = 0.D+00
      ENDDO

       Do J = 1, maxcpus+1
         meandist(J) = 0.00D+00
         sumdist(J) = 0.00D+00
       EndDo

C    find average for each bin across all nodes

C      DO J = 1, NOB - 1
       DO J = 1, NOB
         NDSUM = 0.D+00
         DO I = 1, NNODES
            NDSUM = NDSUM + TESTME(J,I)
         ENDDO
         NDA(J) = NDSUM/DBLE(NNODES)
      ENDDO

C     find sum of distributions
C     find mean of distribution along coordinate
      DO J = 1, NNODES
        sumdist(J) = 0.00D+00
        meandist(J) = 0.00D+00
        Do I = 1, NOB
          sumdist(J) = sumdist(J) + TESTME(I,J)
          meandist(J) = meandist(J) + (TESTME(I,J)*RB(I))
        EndDo
      EndDo
      Do J = 1, NNODES
        meandist(J) = meandist(J)/sumdist(J)
      EndDo
      sumdist(maxcpus+1) = 0.00D+00
      meandist(maxcpus+1) = 0.00D+00
      if(NNODES .eq. 1) then
      Do I = 1, NOB
        sumdist(maxcpus+1) = sumdist(maxcpus+1) +NDA_OLD(I)
        meandist(maxcpus+1) = meandist(maxcpus+1)+(NDA_OLD(I)*RB(I))
      EndDo
      else
      Do I = 1, NOB
        sumdist(maxcpus+1) = sumdist(maxcpus+1) + NDA(I)
        meandist(maxcpus+1) = meandist(maxcpus+1) + (NDA(I)*RB(I))
      EndDo
      endif
      meandist(maxcpus+1) = meandist(maxcpus+1)/sumdist(maxcpus+1)

C    find variance for each bin across all nodes
C      DO J = 1, NOB - 1
      IF(NNODES .eq. 1) then
      DO J = 1, NOB
         NDSIGSUM = 0.D+00
         NDSIGSUM = NDSIGSUM + (NDA(J) - NDA_OLD(J))**2.D+00
         if (NDA(J) .eq. 0.D+00) then
           NDSIG(J) = DSQRT(NDSIGSUM/DBLE(NNODES))
         else
           NDSIG(J) = DSQRT(NDSIGSUM/DBLE(NNODES))/NDA(J)
         endif
      ENDDO
      ELSE
      DO J = 1, NOB
         NDSIGSUM = 0.D+00
C  Change the variance calculation to use nda_old
         DO I = 1, NNODES
            NDSIGSUM = NDSIGSUM + (TESTME(J,I) - NDA(J))**2.D+00
         ENDDO
         if (NDA(J) .eq. 0.D+00) then
           NDSIG(J) = DSQRT(NDSIGSUM/DBLE(NNODES-1))
         else
           NDSIG(J) = DSQRT(NDSIGSUM/DBLE(NNODES-1))/NDA(J)
         endif
      ENDDO
      ENDIF

C     FInd Variance in the distribution means
      NDSIGSUM = 0.00D+00
      Do J = 1, NNODES
       NDSIGSUM = NDSIGSUM + (meandist(J)-meandist(maxcpus+1))**2.D+00
      EndDo
      if (NNODES .gt. 1) then
       NDSIGSUM = NDSIGSUM/DBLE(NNODES-1)
      endif
      NDSIGSUM = DSQRT(NDSIGSUM)/((RB(2) - RB(1))*DBLE(NOB))

CC     Find Maximum Variance
C      NDSGMX = NDSIG(1) 
C      JSTAR = 1 
C
C      DO J = 2, NOB 
C         IF((NDSIG(J).GT.NDSIG(J-1)).AND.(NDSIG(J).GT.NDSGMX))THEN
C            NDSGMX = NDSIG(J)
C            JSTAR = J
C         ENDIF
C      END DO
C
C      IF(NDA(JSTAR) .eq. 0.D+00) then
C      PERCENT = 100.D+00
C      ELSE
C      PERCENT = (NDSGMX)*100.D+00
C      endif

C     Find Percent standard deviation in the means
      PERCENT = NDSIGSUM*100.00D+00

      Do J = 1, NOB
        NDSIG(J) = NDSIG(J)*NDA(J)
      End Do

      RETURN
      END
C***********************************************************************
C  Added NDA_OLD array
C End of SUBROUTINE dntmc_AVERAGEDIST(NNODES,PERCENT,nda_old)
C***********************************************************************

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
C End of double precision FUNCTION dntmc_RAN3(IDUM)
C***********************************************************************

C***********************************************************************
      SUBROUTINE dntmc_NUMBERAVERAGE(L,temp,rconf)
C***********************************************************************
      IMPLICIT none
#include "gibbsP.fh"
#include "util_sgroup.fh"

      integer l   ! input

      integer number
      double precision NUMBERAVR,NUMBERSUM
      double precision rb, testme
      integer  i, j 
      double precision e,k,esqrd,kb,PI
      double precision temp, hold
      logical rconf

      COMMON/NUMBER1/NUMBER(NOB)
      COMMON/NUMBER2/NUMBERSUM(NOB),NUMBERAVR(NOB)
      COMMON/DIST1/RB(NOB),TESTME(NOB,maxcpus)
      COMMON/CONSTANTS/e,k,esqrd,kb,PI

      IF(L.EQ.1)THEN
        DO I = 1, NOB
           NUMBERSUM(I) = 0.D+00
           NUMBERAVR(I) = 0.D+00
        END DO
      ENDIF

      if (rconf) then
        DO I = 1, NOB
          NUMBERSUM(I) = DBLE(NUMBER(I))/DBLE(L)
          NUMBERAVR(I) = NUMBERSUM(I)
          if (i .gt. 1) then
            NUMBERSUM(I) = NUMBERSUM(I) + NUMBERSUM(I-1)
          endif
        enddo
      else
        DO I = 1, NOB
          NUMBERSUM(I) = DBLE(NUMBER(I))/DBLE(L)
          NUMBERAVR(I) = NUMBERSUM(I)
          NUMBERSUM(I) = NUMBERSUM(I)/(4.D+00*PI*(RB(I)**2.D+00))
        END DO
        hold = NUMBERSUM(1)
        DO I = 2, NOB
          if (NUMBERSUM(I) .ge. hold) then
            hold = NUMBERSUM(I)
          endif
        End Do
        DO I = 1, NOB
          NUMBERSUM(I) = NUMBERSUM(I)/hold
        End Do
      endif

      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_NUMBERAVERAGE(L,temp,rconf)
C***********************************************************************

C***********************************************************************
      SUBROUTINE dntmc_NUMBERHIST(L,rconf)
C***********************************************************************
      IMPLICIT none
#include "gibbsP.fh"

      integer l  ! input
      logical rconf ! input

      double precision rmin, rsim, rmax
      integer number

      integer  i
      double precision INCR
      
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/NUMBER1/NUMBER(NOB)

      INCR = (RSIM - RMIN)/(RMAX - RMIN)
      IF(L.EQ.1)THEN
         DO I = 1, NOB
            NUMBER(I) = 0
         END DO
      endif
      
      IF(rconf) then
      I = INT(INCR*DBLE(NOB)) + 1
      else
      I = NINT(INCR*DBLE(NOB))
      endif

      IF (I .le. 0) I = 1
      IF (I .gt. NOB) I = NOB

        NUMBER(I) = NUMBER(I) + 1

      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_NUMBERHIST(L,rconf)
C***********************************************************************

C***********************************************************************
      SUBROUTINE dntmc_MOVE(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
     +               ,ACCEPTPI,UOLD, rtdb, nmols, movetry)
C***********************************************************************
      IMPLICIT none
#include "global.fh"
#include "tcgmsg.fh"
#include "gibbsP.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "util.fh"

      double precision temp  ! input
      LOGICAL ACCEPTT        ! output
      LOGICAL ACCEPTR        ! output
      LOGICAL ACCEPTV        ! output
      LOGICAL ACCEPTPR       ! output
      LOGICAL ACCEPTPI       ! output
      double precision uold  ! input/output
      integer rtdb           ! input
      integer nmols          ! input
      integer movetry        ! output

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
      character*16 tag
      integer ti
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
      COMMON/ATTAG/tag(mxatom*mxmols)      
 
      DIMENSION DELTA(3)

      COMP2 = 0.D+00
      COMP3 = 0.D+00
      COMP2R = 0.D+00
      COMP3R = 0.D+00
      RSTEP = (RMAX - RMIN)
      movetry = 0

c Do these need to more dependant on the actual temperature
c of the simulation?
      PRESS = 1.7D-05   !kcal/mole/A^3  P = 888 torr  taken from UMR data

      RSIMSQ = RSIM*RSIM

C 50% Chance to Change Volume
      if(dntmc_RAN3(IDUM) .lt. 0.5D+00) then

      movetry = movetry + 1
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

      else
      ACCEPTV = .false.
      endif ! end of volume change

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


C 50% Chance to do translation
      IF (dntmc_RAN3(IDUM) .lt. 0.5D+00) then 

      movetry = movetry + 2

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
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 'Energy did not converge'
          endif
          unew = -0.5                        ! set the energy arbitrarily high - H atom
        else
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*)
     & 'Energy converged on the second try'
          endif
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

      if(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_low)) then
      write(12,*) 'coordinates with an energy evaluation'
      endif
      ti = 1
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            if(util_print('MCout',print_low)) then
            WRITE(12,*)tag(ti),(R(I,J,IR),IR=1,3)
            endif
C Added write for configuration files
            if(util_print('ALLdata',print_low)) then
            WRITE(10,*)tag(ti),(R(I,J,IR),IR=1,3)
            endif
            ti = ti + 1
C End Added write
         END DO
      END DO
c
c convert to kcal/mol
c
C Added write for configuration files
      if(util_print('ALLdata',print_low)) then
      WRITE(10,*)'Ucalc = ', Unew,' hartree'
      WRITE(10,*)'Rsim = ', RSIM,' Angstrom'
      call flush(10)
      endif
      endif
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
      if(util_print('ALLdata',print_low)) then
      if(ga_nodeid() .eq. 0) WRITE(10,*)
     &  'Rsim (Calculated) = ',DSQRT(DIST),' Angstrom'
      if(ga_nodeid() .eq. 0) call flush(10)
      endif
C End added write
      unew = unew*627.5095D+00
      if(util_print('MCout',print_low)) then
      if(ga_nodeid() .eq. 0) write(12,*) 'unew = ',unew
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

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
      Else
      ACCEPTT = .false.
      endif ! End of translation
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


C 50% Chance of doing rotation
      if (dntmc_RAN3(IDUM) .lt. 0.5D+00) then

      movetry = movetry + 4

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
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 'Energy did not converge'
          endif
          unew = -0.5                        ! set the energy arbitrarily high - H atom
        else
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*)
     &   'Energy converged on the second try'
          endif
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
      if(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_low)) then
      write(12,*) 'coordinates with an energy evaluation'
      endif
      ti = 1
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            if(util_print('MCout',print_low)) then
            WRITE(12,*)tag(ti),(R(I,J,IR),IR=1,3)
            endif
C Added write for configuration files
            if(util_print('ALLdata',print_low)) then
            WRITE(10,*)tag(ti),(R(I,J,IR),IR=1,3)
            endif
            ti = ti + 1
C End Added write
         END DO
      END DO
c
c convert to kcal/mol
c
C Added write for configuration files
      if(util_print('ALLdata',print_low)) then
      WRITE(10,*)'Ucalc = ', Unew,' hartree'
      WRITE(10,*)'Rsim = ', RSIM,' Angstrom'
      call flush(10)
      endif
      endif
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
      if(util_print('ALLdata',print_low)) then
      if(ga_nodeid() .eq. 0) WRITE(10,*)
     &  'Rsim (Calculated) = ',DSQRT(DIST),' Angstrom'
      if(ga_nodeid() .eq. 0) call flush(10)
      endif
C End added write
      unew = unew*627.5095D+00
      if(util_print('MCout',print_low)) then
      if(ga_nodeid() .eq. 0) write(12,*) 'unew = ',unew
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

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

      ELSE
      ACCEPTR = .false.
      ENDIF ! End of Rotation 

C Added write for configuration files
      if(ga_nodeid() .eq. 0) then
      if(util_print('MCdata',print_low)) then
      ti = 1
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            WRITE(11,*)tag(ti),(R(I,J,IR),IR=1,3)
            ti = ti + 1
         enddo
      enddo
      call flush(11)
      endif
      endif
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
         if(ga_nodeid() .eq. 0) then
         if(util_print('MCdata',print_low)) then
         WRITE(11,*)'Ucalc = ',(Uold/627.5095D+00),' hartree'
         WRITE(11,*)'Rsim = ',rsim,' Angstrom'
         WRITE(11,*)'Rsim (Calculated) = ',DSQRT(DIST),' Angstrom'
         call flush(11)
         endif
         endif
C  End added write
      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_MOVE(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
C     +               ,ACCEPTPI,UOLD, rtdb, nmols, movetry)
C***********************************************************************

C***********************************************************************
      SUBROUTINE dntmc_MOVE_RCONF(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,
     +               ACCEPTPR,ACCEPTPI,UOLD, rtdb, nmols,MOVETRY,dipold)
C***********************************************************************
      IMPLICIT none
#include "global.fh"
#include "tcgmsg.fh"
#include "gibbsP.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "util.fh"

      double precision temp  ! input
      LOGICAL ACCEPTT        ! output
      LOGICAL ACCEPTR        ! output
      LOGICAL ACCEPTV        ! output
      LOGICAL ACCEPTPR       ! output
      LOGICAL ACCEPTPI       ! output
      double precision uold  ! input/output
      integer rtdb           ! input
      integer nmols          ! input
      integer movetry        ! output
      double precision dipold(3) ! input

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
      double precision dipnew(3)
C Added array for configuration files
      double precision DISTcalc(MXMOLS)
      character*16 tag
      integer ti
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
      COMMON/ATTAG/tag(mxatom*mxmols)

      DIMENSION DELTA(3)

      COMP2 = 0.D+00
      COMP3 = 0.D+00
      COMP2R = 0.D+00
      COMP3R = 0.D+00
      movetry = 0

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

C Do translation 50% of the time
      IF (dntmc_RAN3(IDUM) .lt. 0.5D+00) then

      movetry = movetry + 2

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

C      if ((rnew .gt. rmax).or.(rnew .lt. rmin)) goto 20
       if (rnew .gt. rmax) goto 20   

C Do Metropolis on translatonal move

      call dntmc_setcoords(rtdb, r, nmols, natom, spc)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 'Energy did not converge'
          endif
          unew = -0.5                        ! set the energy arbitrarily high - H atom
        else
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*)
     &      'Energy converged on the second try'
          endif
          if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $      call errquit('dntmc_move_rconf: rtdb-task:energy?',0
     &      , RTDB_ERR)
          if (.not. rtdb_get(rtdb, 'task:dipole', mt_dbl, 3, dipnew))
     &       then
             dipnew(1) = 0.00D+00
             dipnew(2) = 0.00D+00
             dipnew(3) = 0.00D+00
          endif

C     &      call errquit('dntmc_move_rconf: rtdb-task:dipole?',0
C     &      , RTDB_ERR)
        endif
      else
        if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $     call errquit('dntmc_move_rconf: rtdb-task:energy?',0
     &     , RTDB_ERR)
          if (.not. rtdb_get(rtdb, 'task:dipole', mt_dbl, 3, dipnew))
     &       then
             dipnew(1) = 0.00D+00
             dipnew(2) = 0.00D+00
             dipnew(3) = 0.00D+00
          endif

c     &      call errquit('dntmc_move_rconf: rtdb-task:dipole?',0
c     &      , RTDB_ERR)
      endif
c      if (.not. task_energy(rtdb)) call errquit
c     $     ('dntmc_move: energy failed',0, GEOM_ERR)
c      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
c     $     call errquit('dntmc_move: rtdb?',0, RTDB_ERR)

      if(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_low)) then
      write(12,*) 'coordinates with an energy evaluation'
      endif
      ti = 1
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            if(util_print('MCout',print_low)) then
            WRITE(12,*)tag(ti),(R(I,J,IR),IR=1,3)
            endif
C Added write for configuration files
            if(util_print('ALLdata',print_low)) then
            WRITE(10,*)tag(ti),(R(I,J,IR),IR=1,3)
            endif
            ti = ti + 1
C End Added write
         END DO
      END DO

c convert to kcal/mol
c
C Added write for configuration files
      if(util_print('ALLdata',print_low)) then
      WRITE(10,*)'Ucalc = ', Unew,' hartree'
      WRITE(10,*)'Dipole = ',(dipnew(IR), IR=1,3),' au'
      WRITE(10,*)'Rsim = ', RNEW,' Angstrom'
      call flush(10)
      endif
      endif
C End added write
      unew = unew*627.5095D+00
      if(util_print('MCout',print_low)) then
      if(ga_nodeid() .eq. 0) write(12,*) 'unew = ',unew
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

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

      IF (TMP.LT.RANN) GOTO 20


C     ACCEPT

  15  CONTINUE

      ACCEPTT = .TRUE.
      UOLD = UNEW
      rsim = rnew
      DO I = 1, 3
        dipold(i) = dipnew(i)
      enddo

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

      else
      ACCEPTT = .false.
      Endif ! End of Translation

C Save Coordinates
      DO I = 1, NMOLS
           DO J = 1, NATOM(SPC(I))
              RS(I,J,1) = R(I,J,1)
              RS(I,J,2) = R(I,J,2)
              RS(I,J,3) = R(I,J,3)
           END DO
      END DO

C Do rotation half the time
      if(dntmc_RAN3(IDUM) .lt. 0.5D+00) then

      movetry = movetry + 4

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

C      if ((rnew .gt. rmax).or.(rnew .lt. rmin)) goto 35
      if (rnew .gt. rmax) goto 35  

C Do Metropolis on rotation

      call dntmc_setcoords(rtdb, r, nmols, natom, spc)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 'Energy did not converge'
          endif
          unew = -0.5                        ! set the energy arbitrarily high - H atom
        else
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*)
     &      'Energy converged on the second try'
          endif
          if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $      call errquit('dntmc_move_rconf: rtdb-task:energy?',0
     &      , RTDB_ERR)
          if (.not. rtdb_get(rtdb, 'task:dipole', mt_dbl, 3, dipnew))
     &       then
             dipnew(1) = 0.00D+00
             dipnew(2) = 0.00D+00
             dipnew(3) = 0.00D+00
          endif

C     &      call errquit('dntmc_move_rconf: rtdb-task:dipole?',0
C     &      , RTDB_ERR)
        endif
      else
        if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
     $     call errquit('dntmc_move_rconf: rtdb-task:energy?',0
     &     , RTDB_ERR)
          if (.not. rtdb_get(rtdb, 'task:dipole', mt_dbl, 3, dipnew))
     &       then
             dipnew(1) = 0.00D+00
             dipnew(2) = 0.00D+00
             dipnew(3) = 0.00D+00
          endif
C     &      call errquit('dntmc_move_rconf: rtdb-task:dipole?',0
C     &      , RTDB_ERR)
      endif
c      if (.not. task_energy(rtdb)) call errquit
c     $     ('dntmc_move: energy failed',0, GEOM_ERR)
c      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, unew))
c     $     call errquit('dntmc_move: rtdb?',0, RTDB_ERR)

      if(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_low)) then
      write(12,*) 'coordinates with an energy evaluation'
      endif
      ti = 1
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            if(util_print('MCout',print_low)) then
            WRITE(12,*)tag(ti),(R(I,J,IR),IR=1,3)
            endif
C Added write for configuration files
            if(util_print('ALLdata',print_low)) then
            WRITE(10,*)tag(ti),(R(I,J,IR),IR=1,3)
            endif
            ti = ti + 1
C End Added write
         END DO
      END DO
c
c convert to kcal/mol
c
C Added write for configuration files
      if(util_print('ALLdata',print_low)) then
      WRITE(10,*)'Ucalc = ', Unew,' hartree'
      WRITE(10,*)'Dipole = ',(dipnew(IR),IR=1,3),' au'
      WRITE(10,*)'Rsim = ', RNEW,' Angstrom'
      call flush(10)
      endif
      endif
C End added write
      unew = unew*627.5095D+00
      if(util_print('MCout',print_low)) then
      if(ga_nodeid() .eq. 0) write(12,*) 'unew = ',unew
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

c      WRITE(12,*)'after rotation from psum'
c      WRITE(12,*)'UNEW =',UNEW
c      WRITE(12,*)'UOLD =',UOLD

      DELTAU = UNEW - UOLD

      ARG = DELTAU/(kb*TEMP)

      TMP = DEXP(-ARG)


      RANN = dntmc_RAN3(IDUM)

      IF (TMP.GE.1.D+00) GOTO 30

      IF (TMP.LT.RANN) GOTO 35

C Accept
   30 CONTINUE

      ACCEPTR = .TRUE.
      UOLD = UNEW
      RSIM = RNEW
      DO I = 1, 3
        dipold(i) = dipnew(i)
      enddo

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

      ELSE
      ACCEPTR = .false.
      ENDIF ! End of Rotation 

C  Finish up
      if(ga_nodeid() .eq. 0) then
C Added write for configuration files
      if(util_print('MCdata',print_low)) then
      ti = 1
      DO I = 1, NMOLS  ! This should be seperated into two loops, copy and print
         DO J = 1, NATOM(SPC(I))
            WRITE(11,*)tag(ti),(R(I,J,IR),IR=1,3)
            ti = ti + 1
         enddo
      enddo
         WRITE(11,*)'Ucalc = ',(Uold/627.5095D+00),' hartree'
         WRITE(11,*)'Dipole =',(dipold(IR),IR=1,3),' au'
         WRITE(11,*)'Rsim = ',rsim,' Angstrom'
         call flush(11)
      endif
      endif
C  End added write
      RETURN
      END
C***********************************************************************
C  End of SUBROUTINE dntmc_MOVE_RCONF(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,
C     +               ACCEPTPR,ACCEPTPI,UOLD, rtdb,nmols,movetry,dipold)
C***********************************************************************

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
C End of SUBROUTINE dntmc_ROOTMEANSQUARE(L,NMOLS,NSPC)
C***********************************************************************

C***********************************************************************
      SUBROUTINE dntmc_AVERAGE(L,TEMP,UU,UAVR,Cv,RSIMAVR,usqrd,rsimsum)
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

      double precision rsimsum, usqrd ! output
      double precision tsqrd

      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/AVR/USUM,USQRDSUM
      COMMON/RADIUS/RMIN,RSIM,RMAX

      IF (L.EQ.1) THEN
         USUM = 0.D+00
         RSIMSUM = 0.D+00
         USQRDSUM = 0.D+00
      ENDIF

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
C End of SUBROUTINE dntmc_AVERAGE(L,TEMP,UU,UAVR,Cv,RSIMAVR,usqrd,rsimsum)
C***********************************************************************

C***********************************************************************
      Subroutine dntmc_read_dntmc_nwchem(rtdb, nmols, prcfil)
C***********************************************************************
      Implicit none
C Include Statements
#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"

C Variable Declarations
C All Internal
      integer ncent, i, j, k, nmols, natom, spc, rdim ! nmols input
      integer rtdb ! input
      double precision r, coords, rcm, rs
      dimension coords(3, MXATOM*MXMOLS)
      character*4 prcfil ! input
      character*16 tag

C Common Blocks
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)
      COMMON/ATTAG/tag(mxatom*mxmols)

C Main Program
      ncent = 0
      Do i = 1, nmols
        Do j = 1, natom(spc(i))
          ncent = ncent + 1
        end do
      end do

      rdim = 3*ncent

      if (.not. rtdb_get(rtdb, 'dntmc:r:'//prcfil(1:4), MT_DBL, rdim,
     &   coords))
     &   call errquit('dntmc_read_dntmc_nwchem:  rtdb error for r',0,
     &   RTDB_ERR)

      if (.not. rtdb_cget(rtdb, 'dntmc:attag:'//prcfil(1:4),
     &   ncent,tag))
     &   call errquit('dntmc_read_dntmc_nwchem: rtdb error for attag',0,
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

C End of subroutine
      return
      end
C***********************************************************************
C  End of Subroutine dntmc_read_dntmc_nwchem(rtdb, nmols, prcfil)
C***********************************************************************

C***********************************************************************
      Subroutine dntmc_readrestart(rtdb, uu, uavr, usqrd, dipole,
     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, ivtry, ittry, irtry)
C***********************************************************************
      Implicit none
C Include Statements
#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"

C Variable Declarations
C All Output
      double precision uu, usum, uavr, usqrdsum, usqrd
      double precision rsim, rsimsum, rsimavr
      double precision seed, r, numbersum, numberavr
      double precision rmax, rmin, rcm, rs
      integer number, iacceptt, iacceptv, iacceptr, iacceptpr, iacceptpi
      integer ivtry, ittry, irtry
      integer rsmcstps
      double precision rsh
      integer rtdb ! input
      integer  nmols !input
      character*4 prcfil !input
      double precision dipole(3) !input

C Common Blocks
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/AVR/USUM,USQRDSUM
      COMMON/NUMBER1/NUMBER(NOB)
      COMMON/NUMBER2/NUMBERSUM(NOB),NUMBERAVR(NOB)
      COMMON/SEEDS/SEED(1)

C Main Program
      if (.not. rtdb_get(rtdb, 'dntmc:seed:'//prcfil(1:4), MT_DBL, 1,
     &   SEED(1)))
     &   call errquit('dntmc_readrestart:  rtdb error for seed(1)',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:uu:'//prcfil(1:4), MT_DBL, 1, UU))
     &   call errquit('dntmc_readrestart:  rtdb error for UU',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:usum:'//prcfil(1:4), MT_DBL, 1,
     &   usum))
     &   call errquit('dntmc_readrestart:  rtdb error for usum',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:uavr:'//prcfil(1:4), MT_DBL, 1,
     &   uavr))
     &   call errquit('dntmc_readrestart:  rtdb error for uavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:usqrdsum:'//prcfil(1:4), MT_DBL, 1
     &   ,usqrdsum))
     &   call errquit('dntmc_readrestart:  rtdb error for usqrdsum',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:usqrd:'//prcfil(1:4), MT_DBL, 1,
     &   usqrd))
     &   call errquit('dntmc_readrestart:  rtdb error for usqrd',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:dipole:'//prcfil(1:4), MT_DBL, 3,
     &   dipole))
     &   call errquit('dntmc_readrestart:  rtdb error for dipole',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rmin:'//prcfil(1:4), MT_DBL, 1,
     &   rmin))
     &   call errquit('dntmc_readrestart:  rtdb error for rmin',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rmax:'//prcfil(1:4), MT_DBL, 1,
     &   rmax))
     &   call errquit('dntmc_readrestart:  rtdb error for rmax',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsim:'//prcfil(1:4), MT_DBL, 1,
     &   rsim))
     &   call errquit('dntmc_readrestart:  rtdb error for rsim',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsimsum:'//prcfil(1:4), MT_DBL, 1,
     &   rsimsum))
     &   call errquit('dntmc_readrestart:  rtdb error for rsimsum',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsimavr:'//prcfil(1:4), MT_DBL, 1,
     &   rsimavr))
     &   call errquit('dntmc_readrestart:  rtdb error for rsimavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:iacceptt:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptt))
     &   call errquit('dntmc_readrestart:  rtdb error for iacceptt',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:iacceptr:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptr))
     &   call errquit('dntmc_readrestart:  rtdb error for iacceptr',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:iacceptv:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptv))
     &   call errquit('dntmc_readrestart:  rtdb error for iacceptv',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:ivtry:'//prcfil(1:4), MT_INT, 1
     &   ,ivtry))
     &   call errquit('dntmc_readrestart:  rtdb error for ivtry',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:ittry:'//prcfil(1:4), MT_INT, 1
     &   ,ittry))
     &   call errquit('dntmc_readrestart:  rtdb error for ittry',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:irtry:'//prcfil(1:4), MT_INT, 1
     &   ,irtry))
     &   call errquit('dntmc_readrestart:  rtdb error for irtry',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:iacceptpr:'//prcfil(1:4), MT_INT,
     &   1, iacceptpr))
     &   call errquit('dntmc_readrestart:  rtdb error for iacceptpr',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:iacceptpi:'//prcfil(1:4), MT_INT,
     &   1, iacceptpi))
     &   call errquit('dntmc_readrestart:  rtdb error for iacceptpi',0,
     &   RTDB_ERR)

      call dntmc_read_dntmc_nwchem(rtdb, nmols, prcfil)

      if (.not. rtdb_get(rtdb, 'dntmc:number:'//prcfil(1:4), MT_INT, NOB
     &   ,number))
     &   call errquit('dntmc_readrestart:  rtdb error for number',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:numbersum:'//prcfil(1:4), MT_DBL,
     &   NOB,
     &   numbersum))
     &   call errquit('dntmc_readrestart:  rtdb error for numbersum',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:numberavr:'//prcfil(1:4), MT_DBL,
     &   NOB,
     &   numberavr))
     &   call errquit('dntmc_readrestart:  rtdb error for numberavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsh:'//prcfil(1:4), MT_DBL, 1, rsh
     &   ))call errquit('dntmc_readrestart:  rtdb error for rsh',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsmcstps:'//prcfil(1:4), MT_INT, 1
     &   ,rsmcstps))
     &   call errquit('dntmc_readrestart:  rtdb error for rsmcstps',0,
     &   RTDB_ERR)

C End Subroutine
      return
      end
C***********************************************************************
C  End of Subroutine dntmc_readrestart(rtdb, uu, uavr, usqrd, dipole,
C     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
C     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, ivtry, ittry, irtry)
C***********************************************************************

C***********************************************************************
      SUBROUTINE dntmc_CMASS(I)
C***********************************************************************
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

      DO J = 1, NATOM(SPC(I))
          RCM(I,1)  = M(SPC(I),J)*R(I,J,1) + RCM(I,1)
          RCM(I,2)  = M(SPC(I),J)*R(I,J,2) + RCM(I,2)
          RCM(I,3)  = M(SPC(I),J)*R(I,J,3) + RCM(I,3)
      END DO

      RCM(I,1) = RCM(I,1) / MTOT(SPC(I))
      RCM(I,2) = RCM(I,2) / MTOT(SPC(I))
      RCM(I,3) = RCM(I,3) / MTOT(SPC(I))

      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_CMASS(I)
C***********************************************************************

C***********************************************************************
      SUBROUTINE dntmc_CLUSTERCMASS(NMOLS)
C***********************************************************************
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

      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_CLUSTERCMASS(NMOLS)
C***********************************************************************

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
C End of SUBROUTINE dntmc_RAN(KK)
C***********************************************************************

C***********************************************************************
      SUBROUTINE dntmc_STRUCT2(K,IMLW,ACCP)
C***********************************************************************
      IMPLICIT none
#include "gibbsP.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"

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

      if(util_print('MCout',print_debug)) then
      if(ga_nodeid() .eq. 0) WRITE(12,*)'d is ', d
      if(ga_nodeid() .eq. 0) call flush(12)
      endif
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

      if(util_print('MCout',print_debug)) then
      if(ga_nodeid() .eq. 0) WRITE(12,*)'structure bad!'
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

      RETURN

   55 CONTINUE

      ACCP = .TRUE.
      
      if(util_print('MCout',print_debug)) then
      if(ga_nodeid() .eq. 0) WRITE(12,*)'structure ok!'
      if(ga_nodeid() .eq. 0) call flush(12)
      endif
 
      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_STRUCT2(K,IMLW,ACCP)
C***********************************************************************

C***********************************************************************
      SUBROUTINE dntmc_MOVEM(K,IMLW,ACCP)
C***********************************************************************
      IMPLICIT none
#include "gibbsP.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"

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
              if(util_print('MCout',print_debug)) then
              if(ga_nodeid() .eq. 0)
     &  WRITE(12,*)'RSAVE(',I,',',J,') =',RSAVE(I,J)
              endif
           END DO
      END DO
      if(util_print('MCout',print_debug)) then
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

C     PICK DISPLACEMENTS ALONG X,Y,Z

      CALL dntmc_RAN(1)
      DLTA(1) = DISPL*(2.D+00*RNNUM(1) - 1.D+00)
      CALL dntmc_RAN(1)
      DLTA(2) = DISPL*(2.D+00*RNNUM(1) - 1.D+00)
      CALL dntmc_RAN(1)
      DLTA(3) = DISPL*(2.D+00*RNNUM(1) - 1.D+00)

      If(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_debug)) then
      DO I = 1, 3
         WRITE(12,*)'DLTA(',I,') =',DLTA(I)
      ENDDO
      call flush(12)
      endif
      endif

C     TRANSLATE COORDINATES OF CHOSEN MOLECULE
c      WRITE(12,*)'coords after translation'

      DO I = 1, NATOM(SPC(K))
         DO J = 1, 3
            RC(K,I,J) = RC(K,I,J) + DLTA(J)
            if(util_print('MCout',print_debug)) then
            if(ga_nodeid() .eq. 0)
     &  WRITE(12,*)'RC(',K,',',I,',',J,') =',RC(K,I,J)
            endif
         END DO
      END DO
      if(ga_nodeid() .eq. 0) call flush(12)


C     PICK AN ANGLE

      CALL dntmc_RAN(1)
      BETA = 2.D+00*ANGLE*(RNNUM(1) - 0.5D+00)

      if(util_print('MCout',print_debug)) then
      if(ga_nodeid() .eq. 0) 
     & WRITE(12,*)'ROTATION ANGLE =',BETA,' radians'
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

      SDBETA = DSIN(BETA)
      CDBETA = DCOS(BETA)

C     PICK AN AXIS ABOUT WHICH TO ROTATE

      CALL dntmc_RAN(1)
      IA1 = INT(3.D+00*RNNUM(1) + 1.D+00)
      IA2 = IA1 + 1
      IA3 = IA1 + 2
      IF(IA2.GT.3) IA2 = IA2 - 3
      IF(IA3.GT.3) IA3 = IA3 - 3

      if(util_print('MCout',print_debug)) then
      if(ga_nodeid() .eq. 0)
     & WRITE(12,*)'ROT.AXIS(xyz->123) =',IA1,IA2,IA3
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

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

      if(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_debug)) then
      DO J = 1, NATOM(SPC(K))
         WRITE(12,*)'RC(',K,',',J,',',1,') =',RC(K,J,1)
         WRITE(12,*)'RC(',K,',',J,',',2,') =',RC(K,J,2)
         WRITE(12,*)'RC(',K,',',J,',',3,') =',RC(K,J,3)
      END DO
      call flush(12)
      endif
      endif

      CALL dntmc_STRUCT2(K,IMLW,ACCP)

c      IF (ACCP) GOTO 12
      IF (ACCP) then
      GOTO 12
      endif

  11  CONTINUE

      ACCP = .FALSE.

      if(util_print('MCout',print_debug)) then
      if(ga_nodeid() .eq. 0) write(12,*) 'inside of 11 continue'
      endif

      DO I = 1, NATOM(SPC(K))
           DO J = 1, 3
               RC(K,I,J) = RSAVE(I,J)
           END DO
      END DO


  12  CONTINUE
      
      if(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_debug)) then
      WRITE(12,*)'after structure routine'
      DO J = 1, NATOM(SPC(K))
         WRITE(12,*)'RC(',K,',',J,',',1,') =',RC(K,J,1)
         WRITE(12,*)'RC(',K,',',J,',',2,') =',RC(K,J,2)
         WRITE(12,*)'RC(',K,',',J,',',3,') =',RC(K,J,3)
      END DO
      call flush(12)
      endif
      endif

      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_MOVEM(K,IMLW,ACCP)
C***********************************************************************

C***********************************************************************
      SUBROUTINE dntmc_CONFIG(IMLW,NODE,RMX,rinit)
C***********************************************************************
      IMPLICIT none
#include "gibbsP.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"

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

      USPC(1) = 2 ! used only for print

      if(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_debug)) then
      WRITE(12,*)'INITIAL CONFIG DATA'
      WRITE(12,*)'IMLW =',IMLW
      WRITE(12,*)'NODE =',NODE
      WRITE(12,*)'RMX =',RMX

      DO I = 1, IMLW
         WRITE(12,*)'SPC(',I,') =',SPC(I)
      ENDDO

      WRITE(12,*)'NATOM(',USPC(1),') =',NATOM(USPC(1)) ! weird print
      call flush(12)
      endif
      endif

      MCSTEPS = 100   !MCSTEPS

      RNODE = DBLE(NODE*NODE)/938.D+00 ! problematic when node > 256

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

c needs to be rinit(spc(i),j,1..3)? Assumes only one species
      DO I = 1, IMLW
         DO J = 1, NATOM(SPC(I))
            RC(I,J,1) = RINIT(1,J,1)
            RC(I,J,2) = RINIT(1,J,2)
            RC(I,J,3) = RINIT(1,J,3)
c            WRITE(12,*)(RC(I,J,IR),IR=1,3)
         END DO
      END DO


      DO K = 1, IMLW

      if(util_print('MCout',print_debug)) then
      if(ga_nodeid() .eq. 0) WRITE(12,*)'MOLECULE =',K
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

      DO LL = 1 , MCSTEPS
         ACCP = .FALSE.

         CALL dntmc_MOVEM(K,IMLW,ACCP)

         IF(ACCP) goto 77
      END DO

  77  CONTINUE

      if(util_print('MCout',print_debug)) then
      if(ga_nodeid() .eq. 0) 
     & WRITE(12,*)'MOLECULE =',K,' WAS MOVED ',LL,' TIMES'
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

      END DO ! K

c      WRITE(12,*)'****************************************'
c      WRITE(12,*)'coordinates from config subroutine...'

      DO I = 1, IMLW
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = RC(I,J,1)
            R(I,J,2) = RC(I,J,2)
            R(I,J,3) = RC(I,J,3)
c           WRITE(12,*)'RC(',I,',',J,',',1,') =',RC(I,J,1)
c           WRITE(12,*)'RC(',I,',',J,',',2,') =',RC(I,J,2)
c           WRITE(12,*)'RC(',I,',',J,',',3,') =',RC(I,J,3)
c           WRITE(12,*)'R(',I,',',J,',',1,') =',R(I,J,1)
c           WRITE(12,*)'R(',I,',',J,',',2,') =',R(I,J,2)
c           WRITE(12,*)'R(',I,',',J,',',3,') =',R(I,J,3)
         END DO
      END DO

c      WRITE(12,*)'****************************************'

      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_CONFIG(IMLW,NODE,RMX,rinit)
C***********************************************************************

C***********************************************************************
      subroutine dntmc_setup_dntmc(rtdb, iacceptt, iacceptr, iacceptv,
     &  ivtry, ittry, irtry,
     &  iacceptpr, iacceptpi, cont1, cont2, nmols,
     &  mprnt, temp, rconf, lmin, lmax, myid, idmedium, groupnums,
     &  prcfil, uu, uavr, usqrd, rsimsum, rsimavr, nspc, dipole,
     &  conv_perc)
C***********************************************************************
      Implicit none

C Include Statements
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "gibbsP.fh"
#include "util.fh"
#include "util_sgroup.fh"

C Variable Declartations
      integer iacceptt, iacceptr, iacceptv, iacceptpr, iacceptpi
      integer ivtry, ittry, irtry
      integer cont1, cont2
      logical os, oc, or, ostate(3) ! internal use only
      integer rtdb ! input
      double precision dipole(3) ! input
      integer nmols
      double precision rinit ! internal use only
      DIMENSION RINIT(MXMOLS,MXATOM,3)
      integer nspc 
      integer mcstps, rsmcstps, mprnt ! mcstps and rsmcstps is internal
      double precision temp
      logical rconf
      integer lmin, lmax
      integer myid, idmedium ! input
      integer node ! internal use only
      integer i, j ! dummy indexes
      double precision rx, ry, rz, dist ! internal use only
      double precision distcalc ! internal use only
      dimension distcalc(MXMOLS) ! internal use only
      integer molmax, IR, natot ! internal use only
      integer groupnums ! input
      dimension groupnums(maxcpus)
      character*4 prcfil ! input
      double precision rsh, rsimsq, rmaxsq ! internal use only
      double precision uu, uavr, usqrd, rsimsum, rsimavr
      logical task_energy
      external task_energy
      external util_get_rtdb_state
      double precision beta, dntmc_RAN3, sdbeta, cdbeta ! internal use only
      double precision comp2, comp3, comp2r, comp3r
      integer ia1, ia2, ia3, ij
      logical new_markov_chain
      double precision radiuscube, rplace, mtotall
      double precision conv_perc

C Variables for Common Blocks
      double precision e, k, esqrd, kb, pi
      double precision disp, alpha
      INTEGER ROTF
      integer idum
      double precision seed
      INTEGER AB
      integer natom
      integer uspc
      double precision m, mtot
      double precision rmin, rmax, rsim
      integer spc
      double precision rclcm
      double precision r, rcm, rs, ri
      integer nmin, ncurr, nmax
      character*16 tag

C Common Blocks
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/TRPS/DISP,ALPHA,ROTF
      COMMON/SEED/IDUM
      COMMON/SEEDS/SEED(1)
      COMMON/AORB/AB
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/USPECIES/USPC(MXSP)
      COMMON/MASSES/M(MXSP,MXATOM),MTOT(MXSP)
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/SPECIES/SPC(MXMOLS)
      COMMON/COORD2/RCLCM(3)
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/COORDI/RI(MXMOLS,MXATOM,3)
      COMMON/GRAND/NMIN,NCURR,NMAX
      COMMON/ATTAG/tag(mxatom*mxmols)


C Main Program
C  Setup Constants
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
      IVTRY = 0
      ITTRY = 0
      IRTRY = 0
C Some constants
      ROTF = 0
      AB = 1
      IDUM = -1
      do i=1, util_sgroup_numgroups()
C        if (groupnums(i) .eq. idmedium) then
        if (groupnums(i) .eq. util_sgroup_mygroup()) then
          IDUM = i
        endif
      enddo
      SEED(1) = 962525439.0

C Aquire Input Information
      if (.not. rtdb_get(rtdb, 'dntmc:norestart', MT_INT, 1, CONT1))
     &   CONT1 = 1

      call util_get_rtdb_state(rtdb, os, oc, or)
      if (or .eq. .true.) then
         CONT2=1
         ostate(1)= os
         ostate(2)= oc
         ostate(3)= .false.
         if(.not. rtdb_put(rtdb,'rtdb:stored:state',mt_log,3,ostate))
     &   call errquit('dntmc_setup_dntmc: rtdb error for restart reset',
     &   0, RTDB_ERR)
      else
         CONT2=0
      endif

      call dntmc_getspecinfo(rtdb, nmols, natom, rinit, uspc, nspc,
     &      m, SPC, tag)

      if (.not.rtdb_get(rtdb, 'dntmc:mcsteps', MT_INT, 1, MCSTPS))
     &    call errquit('dntmc_setup_dntmc: rtdb error for mcsteps', 0,
     &    RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'dntmc:rsmcstps:'//prcfil(1:4),
     &   MT_INT, 1, rsmcstps))
     &   rsmcstps = 0

      if(.not. rtdb_get(rtdb, 'dntmc:mprnt', MT_INT, 1, mprnt))
     &  mprnt = 10
      if (.not. rtdb_get(rtdb, 'dntmc:disp', MT_DBL, 1, disp))
     &   disp = 0.04D+00
      if (.not. rtdb_get(rtdb, 'dntmc:alpha', MT_DBL, 1, alpha))
     &   alpha = 0.06D+00
      if (.not.rtdb_get(rtdb, 'dntmc:temp', MT_DBL, 1, TEMP))
     &    call errquit('dntmc_setup_dntmc: rtdb error for temp',
     &     0, RTDB_ERR)
      if (.not.rtdb_get(rtdb, 'dntmc:rmin', MT_DBL, 1, RMIN))
     &    call errquit('dntmc_setup_dntmc: rtdb error for rmin',
     &     0, RTDB_ERR)
      if (.not.rtdb_get(rtdb, 'dntmc:rmax', MT_DBL, 1, RMAX))
     &    call errquit('dntmc_setup_dntmc: rtdb error for rmax',
     &     0, RTDB_ERR)


      if (.not. rtdb_get(rtdb, 'dntmc:rconf', MT_LOG, 1, rconf))
     &   rconf = .true.

      if (.not. rtdb_get(rtdb, 'dntmc:convergence', MT_DBL,1,conv_perc))
     &   conv_perc = 0.00D+00

C Do some computations
      IF(CONT2.EQ.0)THEN
         LMIN = 1
         LMAX = MCSTPS
      ELSE
         LMIN = rsmcstps
         LMAX = rsmcstps + MCSTPS
      ENDIF
    
      RSIM = RMAX

C Assuming that all species are the same
C      DO I=1, NMOLS
C        SPC(I) = USPC(1)
C      ENDDO

C Print some information
C      IF(ga_nodeid() .eq. 0) then
C      if(util_print('MCout',print_debug)) then
C      WRITE(12,*)NSPC,NMOLS,MCSTPS,MPRNT,TEMP,DISP,ALPHA,ROTF
C      WRITE(12,*)AB
C      WRITE(12,*)'CONT1 =',CONT1
C      WRITE(12,*)'CONT2 =',CONT2
C      call flush(12)
C      endif
C
C      IF(myid .eq. 0) then
C      if (util_print('Information',print_none)) then
C      write(6,*)'NSPC  =',NSPC
C      write(6,*)'NMOLS =',NMOLS
C      write(6,*)'TEMP  =',TEMP,'K'
C      write(6,*)'MCSTPS=',MCSTPS
C      write(6,*)'Trans Max Disp =',DISP,'Ang'
C      write(6,*)'Rot Max Disp =',ALPHA,'rad'
C      if (CONT1 .eq. 0) write(6,*)'No restart files written'
C      if (CONT1 .eq. 1) write(6,*)'Restart files written'
C      if (CONT2 .eq. 1) write(6,*)'Restart Run'
C      if (CONT2 .eq. 0) write(6,*)'Start Run'
C      write(6,*)'RMIN =',RMIN
C      write(6,*)'RMAX =',RMAX
C      if (rconf) write(6,*)'R-Config Methodology'
C      if (.not. rconf) write(6,*)'R-Sim Methodology'
C      write(6,*)'NGROUPS =',util_sgroup_numgroups()
C      WRITE(6,*)NSPC,NMOLS,MCSTPS,MPRNT,TEMP,DISP,ALPHA,ROTF
C      WRITE(6,*)AB
C      WRITE(6,*)'CONT1 =',CONT1
C      WRITE(6,*)'CONT2 =',CONT2
C      call flush(6)
C      endif
C      endif
C      endif

C Mass of the system
      DO I = 1, NSPC
         DO J = 1, NATOM(USPC(I))
            MTOT(USPC(I)) = MTOT(USPC(I)) + M(USPC(I),J)
      IF(ga_nodeid() .eq. 0) then
            if(util_print('MCout',print_debug)) then
            WRITE(12,*)'M(',USPC(I),',',J,') =',M(USPC(I),J)
            call flush(12)
            endif
      endif
         END DO
      END DO

      MTOTALL = 0.00D+00
      DO I = 1, NMOLS
        MTOTALL = MTOTALL + MTOT(SPC(I))
      ENDDO

C Generate Initial Configuration if 1st run
      IF(CONT2 .eq. 0) then

      NODE = -1
      do i=1, util_sgroup_numgroups()
C        if (groupnums(i) .eq. idmedium) then
        if (groupnums(i) .eq. util_sgroup_mygroup()) then

          NODE = i
        endif
      enddo
      radiuscube = RMIN
      radiuscube = (DBLE(NODE)*(RMAX-RMIN)/DBLE(util_sgroup_numgroups()
     & +1)) + radiuscube
      rplace = radiuscube
C      CALL dntmc_CONFIG(NMOLS,NODE,RMAX,rinit)
C      CALL dntmc_initialstate(nmols, rinit)
      CALL dntmc_initialstate2(nmols, rinit, radiuscube)
C      NODE = idmedium
       NODE = util_sgroup_mygroup()
C Center on COM
      CALL dntmc_CLUSTERCMASS(NMOLS)
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
      do i=1, util_sgroup_numgroups()
C        if (groupnums(i) .eq. idmedium) then
        if (groupnums(i) .eq. util_sgroup_mygroup()) then

          NODE = i
        endif
      enddo
      DIST = DSQRT(DIST)
      RX = RCM(molmax,1)/DIST
      RY = RCM(molmax,2)/DIST
      RZ = RCM(molmax,3)/DIST

      if (DIST .gt. radiuscube) goto 99
      if (.not. DIST .gt. rplace) then
C Distribute Evenly 100%
      DIST = RMIN - DIST
      DIST = (DBLE(NODE)*(RMAX-RMIN)/DBLE(util_sgroup_numgroups()+1)) 
     &       + DIST
      DIST = DIST / (1.00D+00 - (MTOT(SPC(molmax))/MTOTALL))
C      NODE = idmedium
      NODE = util_sgroup_mygroup()

C Translate by constant
      DO I = 1, NATOM(SPC(molmax))
            R(molmax,I,1) = R(molmax,I,1)+(RX * DIST)
            R(molmax,I,2) = R(molmax,I,2)+(RY * DIST)
            R(molmax,I,3) = R(molmax,I,3)+(RZ * DIST)
      END DO
      endif

CC Do twice
C      Do ij = 1, 2
CC Randomly Rotate Monomers
C      DO I = 1, NMOLS
CC Angle 
C        Beta = 3.14D+00*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
C        SDBETA = DSIN(Beta)
C        CDBETA = DCOS(Beta)
CC Axis
C        IA1 = INT(3.D+00*dntmc_RAN3(IDUM) + 1)
C        IA2 = IA1 + 1
C        IA3 = IA1 + 2
C        IF(IA2.GT.3) IA2 = IA2 - 3
C        IF(IA3.GT.3) IA3 = IA3 - 3
CC Rotate
C        call dntmc_CMASS(I)
C        Do J = 1, NATOM(SPC(I))
C            COMP2 = (R(I,J,IA2) - RCM(I,IA2))
C            COMP3 = (R(I,J,IA3) - RCM(I,IA3))
C
C            COMP2R=CDBETA*COMP2+SDBETA*COMP3
C            COMP3R=-SDBETA*COMP2+CDBETA*COMP3
C
C            R(I,J,IA2)= RCM(I,IA2) + COMP2R
C            R(I,J,IA3)= RCM(I,IA3) + COMP3R
C        END DO
C      END DO
C      END DO
C        
C Recenter
      CALL dntmc_CLUSTERCMASS(NMOLS)
      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO
C set rsim if rconf method
      if(rconf) then
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

C Write Varification
      if(ga_nodeid().eq. 0) then
      if(util_print('MCout',print_debug)) then
      WRITE(12,*)'made it through the initial config setup'
      WRITE(12,*)'on group =',prcfil(1:4)
      call flush(12)
      endif
      if (myid .eq. 0) then
      if(util_print('Information',print_none)) then
      WRITE(6,'(a)')'Made it through the initial config setup'
      call flush(6)
      endif
      endif
      endif

C Else a restart run
      ELSE
        if(ga_nodeid() .eq. 0) then
        if(util_print('Information',print_none)) then
        if(myid .eq. 0) then
        write(6,'(a)')'Reading Restart Information'
        call flush(6)
        endif
        endif
        if(util_print('MCout',print_debug)) then
        write(12,*)'Reading Restart Information for group ',prcfil(1:4)
        call flush(12)
        endif
        endif
        call dntmc_readrestart(rtdb, uu, uavr, usqrd, dipole,
     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, ivtry, ittry, irtry)
        if(ga_nodeid() .eq. 0) then
        if(util_print('Information',print_none)) then
        if(myid .eq. 0) then
        write(6,'(a)')'Done with restart setup'
        call flush(6)
        endif
        endif
        if(util_print('MCout',print_debug)) then
        write(12,*)'Done with restart setup for group ',prcfil(1:4)
        call flush(12)
        endif
        endif
CC Redistribute Starting Structures if rtdb flag is set
C      if (.not. rtdb_get(rtdb, 'dntmc:new_markov_chain', MT_log, 1,
C     & NEW_MARKOV_CHAIN)) then
C         NEW_MARKOV_CHAIN = .false.
C      else
C         NEW_MARKOV_CHAIN = .true.
C            if(.not. rtdb_delete(rtdb, 'dntmc:new_markov_chain'))
C     &        call errquit('dntmc_setup_dntmc:rtdb failed to delete
C     &        dntmc:new_markov_chain'
C     &        ,0,RTDB_ERR)
C      endif
C        if(new_markov_chain) then
CC find outermost molecule and vector
C      DO I = 1, NMOLS
C         CALL dntmc_CMASS(I)
C         RX = RCM(I,1)
C         RY = RCM(I,2)
C         RZ = RCM(I,3)
C         DISTcalc(I) = RX*RX + RY*RY + RZ*RZ
C      END DO
C      DIST = DISTcalc(1)
C      molmax = 1
C      DO I = 2, NMOLS
C         IF (DISTcalc(I) .gt. DIST) then
C            DIST = DISTcalc(I)
C            molmax = I
C         endif
C      enddo
C      NODE = -1
C      do i=1, util_sgroup_numgroups()
C        if (groupnums(i) .eq. idmedium) then
C          NODE = i
C        endif
C      enddo
C      DIST = DSQRT(DIST)
C      RX = RCM(molmax,1)/DIST
C      RY = RCM(molmax,2)/DIST
C      RZ = RCM(molmax,3)/DIST
CC Distribute Evenly 100%
C      DIST = RMIN - DIST
C      DIST = (DBLE(NODE)*(RMAX-RMIN)/DBLE(util_sgroup_numgroups()+1)) + DIST
CC      DIST = ((RMAX+RMIN)/2.D+00) - DIST
CC      DIST = (DBLE(NODE)*(RMAX-RMIN)/2.D+00/DBLE(util_sgroup_numgroups()+1)) + DIST
C      if (NMOLS .eq. 1) then
C        DIST = 0.00D+00
C      else
C        DIST = (DBLE(NMOLS)/DBLE(NMOLS -1)) * DIST
C      endif
C      NODE = idmedium
CC Translate by constant
C      DO I = 1, NATOM(SPC(molmax))
C            R(molmax,I,1) = R(molmax,I,1)+(RX * DIST)
C            R(molmax,I,2) = R(molmax,I,2)+(RY * DIST)
C            R(molmax,I,3) = R(molmax,I,3)+(RZ * DIST)
C      END DO
C
CC Do twice
C      Do ij = 1, 2
CC Randomly Rotate Monomers
C      DO I = 1, NMOLS
CC Angle
C        Beta = 3.14D+00*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
C        SDBETA = DSIN(Beta)
C        CDBETA = DCOS(Beta)
CC Axis
C        IA1 = INT(3.D+00*dntmc_RAN3(IDUM) + 1)
C        IA2 = IA1 + 1
C        IA3 = IA1 + 2
C        IF(IA2.GT.3) IA2 = IA2 - 3
C        IF(IA3.GT.3) IA3 = IA3 - 3
CC Rotate
C        call dntmc_CMASS(I)
C        Do J = 1, NATOM(SPC(I))
C            COMP2 = (R(I,J,IA2) - RCM(I,IA2))
C            COMP3 = (R(I,J,IA3) - RCM(I,IA3))
C
C            COMP2R=CDBETA*COMP2+SDBETA*COMP3
C            COMP3R=-SDBETA*COMP2+CDBETA*COMP3
C
C            R(I,J,IA2)= RCM(I,IA2) + COMP2R
C            R(I,J,IA3)= RCM(I,IA3) + COMP3R
C        END DO
C      END DO
C      END DO
C
CC Recenter
C      CALL dntmc_CLUSTERCMASS(NMOLS)
C      DO I = 1, NMOLS
C         DO J = 1, NATOM(SPC(I))
C            R(I,J,1) = R(I,J,1) - RCLCM(1)
C            R(I,J,2) = R(I,J,2) - RCLCM(2)
C            R(I,J,3) = R(I,J,3) - RCLCM(3)
C         END DO
C      END DO
CC set rsim if rconf method
C      if(rconf) then
C      DO I = 1, NMOLS
C         CALL dntmc_CMASS(I)
C         RX = RCM(I,1)
C         RY = RCM(I,2)
C         RZ = RCM(I,3)
C         DISTcalc(I) = RX*RX + RY*RY + RZ*RZ
C      END DO
C      DIST = DISTcalc(1)
C      DO I = 2, NMOLS
C         IF (DISTcalc(I) .gt. DIST) then
C            DIST = DISTcalc(I)
C         endif
C      enddo
C      rsim = DSQRT(DIST)
C      else
C      rsim = rmax
C      endif
C        endif
      ENDIF

C Setup some more constants
      NMIN = NMOLS - 4
      NCURR = NMOLS
      NMAX = NMOLS
      RSIMSQ = RSIM*RSIM
      RMAXSQ = RMAX*RMAX

C Write Initial Conditions
      IF(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_debug)) then
      WRITE(12,*)'RMAXSQ =',RMAXSQ
      WRITE(12,*)'RSIMSQ =',RSIMSQ
      WRITE(12,*)'TEMPERATURE =',TEMP,' K'
      WRITE(12,*)'RMIN =',RMIN,' ANGSTROMS'
      WRITE(12,*)'RMAX =',RMAX,' ANGSTROMS'
      WRITE(12,*)'The initial simulation radius =',RSIM,' ANGSTROMS'
      WRITE(12,*)'NCURR =',NCURR
      call flush(12)
      endif
      IF(myid .eq. 0) then
      if(util_print('Information',print_none)) then
      write(6,'(a,i10)')'NSPC  =',NSPC
      write(6,'(a,i10)')'NMOLS =',NMOLS
      write(6,'(a,F10.2,a)')'TEMP  =',TEMP,' K'
      write(6,'(a,ES10.2E2)')'MCSTPS=',DBLE(MCSTPS)
      write(6,'(a,ES10.2E2,a)')'Trans Max Disp =',DISP,' Ang'
      write(6,'(a,ES10.2E2,a)')'Rot Max Disp =',ALPHA,' rad'
      if (CONT1 .eq. 0) write(6,'(a)')'No restart files written'
      if (CONT1 .eq. 1) write(6,'(a)')'Restart files written'
      if (CONT2 .eq. 1) write(6,'(a)')'Restart Run'
      if (CONT2 .eq. 0) write(6,'(a)')'Start Run'
      write(6,'(a,F10.2)')'RMIN =',RMIN
      write(6,'(a,F10.2)')'RMAX =',RMAX
      if (rconf) write(6,'(a)')'R-Config Methodology'
      if (.not. rconf) write(6,'(a)')'R-Sim Methodology'
      write(6,'(a,i10)')'NGROUPS =',util_sgroup_numgroups()
C      WRITE(6,*)NSPC,NMOLS,MCSTPS,MPRNT,TEMP,DISP,ALPHA,ROTF
C      WRITE(6,*)AB
C      WRITE(6,*)'CONT1 =',CONT1
C      WRITE(6,*)'CONT2 =',CONT2
      call flush(6)
      endif
      ENDIF
      ENDIF

C Check Validity of initial coordinates
C Recenter Cluster
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
         IF(rconf) then
           IF(DIST.GT.RMAXSQ) goto 99
         ELSE
           IF(DIST.GT.RSIMSQ) goto 99
         ENDIF
      END DO
      
      Goto 101

  99  CONTINUE
      IF(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_debug)) then
      WRITE(12,*)'dist is ',DSQRT(DIST),' mol ',i,' atom ',j
C should never get here unless the algorithm is messed up
      WRITE(12,*)'molecule',I,' is out of the simulation volume'
      call flush(12)
      endif
      endif
      call errquit('dntmc_setup_dntmc: should not get to 99', 1,
     &   UNKNOWN_ERR)

 101  CONTINUE

C Save Inital Coordinates
      IF(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_low)) then
      WRITE(12,*)'The initial cluster coordinates for a fixed center'
      WRITE(12,*)'of mass at the origin are:'
      call flush(12)
      endif
      endif

      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            if(util_print('MCout',print_low)) then
            IF(ga_nodeid() .eq. 0) WRITE(12,*)(R(I,J,IR),IR=1,3)
            endif
            RI(I,J,1) = R(I,J,1)
            RI(I,J,2) = R(I,J,2)
            RI(I,J,3) = R(I,J,3)
         END DO
      END DO

      IF(ga_nodeid() .eq. 0) call flush(12)
      call dntmc_setcoords(rtdb, r, nmols, natom, spc)

C Get initial energy if first run or new markov chain
      IF((CONT2 .eq. 0).or.(new_markov_chain)) then
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 'Energy did not converge' 
          endif
          uu = -0.5                        ! set the energy arbitrarily high - H atom
        else
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 
     &    'Energy converged on the second try'
          endif
          if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, uu))
     $      call errquit('dntmc_setup_dntmc: rtdb-task:energy?',0,
     &       RTDB_ERR)
          if (.not. rtdb_get(rtdb, 'task:dipole', mt_dbl, 3, dipole))
     &       then
             dipole(1) = 0.00D+00
             dipole(2) = 0.00D+00
             dipole(3) = 0.00D+00
          endif
C     &      call errquit('dntmc_setup_dntmc: rtdb-task:dipole?',0
C     &      , RTDB_ERR)
        endif
      else
      if (.not. rtdb_get(rtdb, 'task:energy', mt_dbl, 1, uu))
     $     call errquit('dntmc_setup_dntmc: rtdb-task:energy?',0,
     &      RTDB_ERR)
          if (.not. rtdb_get(rtdb, 'task:dipole', mt_dbl, 3, dipole))
     &       then
             dipole(1) = 0.00D+00
             dipole(2) = 0.00D+00
             dipole(3) = 0.00D+00
          endif
C     &      call errquit('dntmc_setup_dntmc: rtdb-task:dipole?',0
C     &      , RTDB_ERR)
      endif

      uu = uu*627.5095D+00
      endif
      if(util_print('MCout',print_low)) then
      IF(ga_nodeid() .eq. 0) WRITE(12,*)'INITIAL UU = ',UU,' kcal/mole'
      if(ga_nodeid() .eq. 0) call flush(12)
      endif

C Number of atoms
      NATOT = 0
      DO I = 1, NMOLS
         NATOT = NATOT + NATOM(SPC(I))
      END DO
C Pring Initialization Varification
      IF(ga_nodeid() .eq. 0) then
      if(util_print('MCout',print_debug)) then
      WRITE(12,*)'LMIN =',LMIN,'  LMAX =',LMAX
      call flush(12)
      endif
      IF(myid .eq. 0) then    
      if(util_print('Information',print_none)) then
      WRITE(6,*)'Made it through initialization'
      call flush(6)
      endif
      ENDIF
      ENDIF
C All Job information is now setup
C End of subroutine
      return
      end

C***********************************************************************
C End of subroutine dntmc_setup_dntmc(rtdb, iacceptt, iacceptr, iacceptv,
C     &  ivtry, ittry, irtry,
C     &  iacceptpr, iacceptpi, cont1, cont2, nmols,
C     &  mprnt, temp, rconf, lmin, lmax, myid, idmedium, groupnums,
C     &  prcfil, uu, uavr, usqrd, rsimsum, rsimavr, nspc, dipole,
C     &  conv_perc)
C***********************************************************************

C***********************************************************************
      subroutine dntmc_setup_group
     &(rtdb, idmedium, myid, groupnums, prcfil, prcfil2, fprefix)
C***********************************************************************
      Implicit none

C Include Statements
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "gibbsP.fh"
#include "util.fh"
#include "util_sgroup.fh"

C Variable declarations
      integer rtdb ! input
      integer myid ! input
      integer idmedium ! input
      integer groupnums ! input
      character*4 prcfil, prcfil2
      integer i
C Indect output of file units
      character*256 fprefix, fprefixcat ! internal use only !fprefix output
      character*256 fprefix2, dntmc_dir

C Dimensions
      dimension groupnums(maxcpus)

C Begin Main Program
      call dntmc_write_prcfil(idmedium, groupnums, prcfil)
      call dntmc_build_prcfil(prcfil2, myid)
      write(fprefix2,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix2))
     &     call errquit('dntmc_setup_group:rtdb get file_prefix failed'
     &     ,0,RTDB_ERR)
      write(dntmc_dir,'(256(a))') (' ', i=1,256)
      if (.not. rtdb_cget(rtdb, 'dntmc:directory',1,dntmc_dir))
     &   dntmc_dir(1:2)='./'
      write(fprefix,'(256(a))') (' ', i=1, 256)
      write(fprefix,'(3(a))')
     &     dntmc_dir(1:inp_strlen(dntmc_dir)),
     &     '/',
     &     fprefix2(1:inp_strlen(fprefix2))

C Open file Units
C general output ! Only group Zeros write
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(3(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.MCout.',
     &      prcfil(1:4)

      if(util_print('MCout',print_low)) then
      OPEN(UNIT = 12,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')
      endif

C  Added files for configurations ! Only group Zeros write
C  Unit 10 ALL
C  Unit 11 Accepted Only
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(3(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.ALLdata.',
     &      prcfil(1:4)
      if(util_print('ALLdata',print_low)) then
      OPEN(UNIT = 10,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')
      endif

      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(3(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.MCdata.',
     &      prcfil(1:4)

      if(util_print('MCdata',print_low)) then
      OPEN(UNIT = 11,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     &    ,STATUS = 'REPLACE')
      endif

C probability distribution function, permanent ! Only node zero writes
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(2(a))')
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.fdist'

      if (myid .eq. 0) then
      if (util_print('fdist',print_low)) then
      OPEN(UNIT = 3,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'REPLACE')
      endif
      endif

C End subroutine
      return
      end

C***********************************************************************
C End subroutine dntmc_setup_group
C     &(rtdb, idmedium, myid, groupnums, prcfil, prcfil2,fprefix)
C***********************************************************************

C***********************************************************************
      subroutine dntmc_setup_ga
     &(rtdb, myid, idbig, inodesbig,
     &idmedium, idzero, g_pr, groupnums)
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
#include "gibbsP.fh"
#include "util.fh"
#include "util_sgroup.fh"

C Variable Declarations
      integer rtdb ! input
      integer myid, idbig, inodesbig ! input
C Depreciated
C      logical util_sgstart
      logical status ! internal use
C Depreciated
C      external util_sgstart
      integer idmedium
      integer ld(2)! internal use
      integer g_proc, procnums, groupnums, bigproc !p_proc is internal
C     procnums is internal
      integer i, j, idzero ! i and j are internal indexes
      integer ndim, dims(2), chunk(2), g_pr !n(dim)s and chunck are internal
C Replacement for util_sgstart
      external util_sggo
C Additional variables
      integer groups_want, array_cpu(1), method
      integer dir
C End Additions

C Dimensions
      dimension procnums(maxcpus)
      dimension groupnums(maxcpus)
      dimension bigproc(maxcpus,1)

C Begin main program
C Create big GAs
      if (.not.ga_create(MT_DBL, inodesbig, 1, "proc list",
     &     1, -1, g_proc))
     &   call errquit('dntmc_setup_ga:g_proc create error', 0,
     &                 GA_ERR)
      call ga_fill(g_proc, -1)

C Create Subgroups
C Depreciated
C      if (.not.util_sgstart(rtdb))
C     &  call errquit('dntmc_setup_ga:error starting subgroups',
C     &                  0, UNKNOWN_ERR)
C Replacement for util_sgstart
      if (.not.rtdb_get(rtdb, 'subgroups_number', mt_int, 1,
     &   groups_want))groups_want=1
C Setting Method = 1
c     1 -- use groups_want to generate equal sized groups (array_cpu ignored)
C     Simplest. 
c     2 -- turn each SMP box into a group (array_cpu and groups_want ignored)
C     This uses GA to tell it about the cluster.
c     3 -- use array_cpu(groups_want) to define number of nodes per group
c     4 -- use array_cpu(groups_want+nnodes) to define which nodes per group
C     This is just option 3, but you get to lay the groups out exactly.
      method = 1
C Setting array_cpu(1) = 0
C Ignored in methods 1 and 2, otherwise must change
      array_cpu(1) = 0
      dir = 1 ! Write group rtdb's in scratch directories
C Replacement call for util_sgstart
      call util_sggo(rtdb,groups_want,method,array_cpu,dir)
C End of changes for util_sgstart depreciation

      if (util_print('debug',print_debug)) then
        write(6,*)'Now in Subgroups'
        call flush(6)
      endif
      idmedium = ga_pgroup_get_default() ! Now in Subgroups
C New version also sets up zero group (util_sgstart)
      idzero = util_sgroup_zero_group()

C Get Back to Big Group
      call ga_pgroup_set_default(idbig) ! Now in Big Group

C Create Processor Zero and Group Processor Lists
      if (ga_pgroup_nodeid(idmedium).eq.0) then
        ld(1)  = 1   !must be physical dimension of local array
        ld(2)  = 1
C        call ga_put(g_proc,myid+1,myid+1,1,1,idmedium,ld)
        call ga_put(g_proc,myid+1,myid+1,1,1,util_sgroup_mygroup(),ld)
        call ga_sync()
        ld(1)  = maxcpus  !must be physical dimension of local array
        ld(2)  = 1
        call ga_get(g_proc,1,inodesbig,1,1,bigproc,ld)
        j = 0
        do i = 1, inodesbig
          if (bigproc(i,1).ne.-1) then
            j = j + 1
            groupnums(j) = bigproc(i,1)
            procnums(j) = i-1
          endif
        enddo
        if (j.ne.util_sgroup_numgroups())
     &    call errquit('dntmc_setup_ga:zero node creation problem',
     &      j, UNKNOWN_ERR)
      else
        call ga_sync()
      endif

C Broadcase Results
      call ga_pgroup_brdcst(idmedium,msg_dntmc3+MSGINT, procnums,
     &     mitob(maxcpus), 0)
      call ga_pgroup_brdcst(idmedium,msg_dntmc7+MSGINT, groupnums,
     &     mitob(maxcpus), 0)

C Destroy GA "proc list"
      if (.not. ga_destroy(g_proc)) call errquit('dntmc_setup_ga:
     &ga_destroy(g_proc) failed', GA_ERR)

C Is done by the new version of util_sgstart
C Create Zero Group
C      idzero = ga_pgroup_create(procnums, util_sgroup_numgroups())

C Create Zeros GAs
      if (ga_pgroup_nodeid(idmedium).eq.0) then
        call ga_pgroup_set_default(idzero)
          ndim = 2
          dims(1) = NOB
          dims(2) = util_sgroup_numgroups()
          chunk(1) = dims(1)
          chunk(2) = 0
          g_pr = ga_create_handle()
          call ga_set_data(g_pr,ndim,dims,MT_F_DBL)
          call ga_set_chunk(g_pr,chunk)
          status = ga_allocate(g_pr)
          call ga_zero(g_pr)
       endif

C Start Subgroups
       call ga_pgroup_set_default(idmedium) ! Now in Subgroups
       call ga_sync()
C Test write
C       write(6,*) '@ I am proc ',myid,' in group ',idmedium
C     & ,' with group proc ',ga_nodeid()
C       write(6,*) '@ I am proc ',myid,' in group2 ',my_ga_grp(depth),
C     & ' with zero group ',util_sgroup_zero_group()
C       if (myid .eq. 0) then
C           write(6,*)'@ group and proc arrays'
C         Do i=1,inodesbig
C           write(6,*)'@ ',groupnums(i),procnums(i)
C         EndDo
C       endif

C End Routine
       return
       end

C***********************************************************************
C End subroutine dntmc_setup_ga
C     &(rtdb, myid, idbig, inodesbig, 
C     &idmedium, idzero, g_pr,  groupnums)
C***********************************************************************

C***********************************************************************
      subroutine dntmc_write_prcfil(idmedium, groupnums, prcfil)
C***********************************************************************
      Implicit none

C Include Statements
#include "gibbsP.fh"
#include "errquit.fh"
#include "util_sgroup.fh"

C Variable Declarations
      integer i
      integer idmedium, groupnums
      integer groupid
      character*4 prcfil ! only output
      dimension groupnums(maxcpus)

C Main Program
      groupid = -1
      do i=1, util_sgroup_numgroups()
C        if (groupnums(i) .eq. idmedium) then
        if (groupnums(i) .eq. util_sgroup_mygroup()) then

          groupid = i
        endif
      enddo

      if (groupid .eq. -1)
     & call errquit('dntmc_write_prcfil:failed to allocate idgroup',0,

     &   GA_ERR)

       call dntmc_build_prcfil(prcfil, groupid)

C End Main Program
      return
      end

C***********************************************************************
C  End of subroutine dntmc_write_prcfil(idmedium, groupnums, prcfil)
C***********************************************************************

C***********************************************************************
      subroutine dntmc_restart_post(procrestart,rtdb)
C***********************************************************************
      Implicit None

C Include Statements
#include "tcgmsg.fh"
#include "global.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"

C Variable Declarations
      integer rtdb, procrestart ! input
      character*256 fprefix
      integer nmols, natom
      double precision rinit(MXMOLS,MXATOM,3)
      integer uspc(MXSP), spc(MXMOLS), nspc
      double precision m(MXMOLS,MXATOM)
      integer i !input
      character*4 prcfil
      character*16 tag(mxatom*mxmols) ! internal
      character*256 dntmc_dir
C Everything else is internal

C Main Program
      write(fprefix,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix))
     &  call errquit('dntmc_restart_post:rtdb get file_prefix failed',0
     &  ,RTDB_ERR)
      write(dntmc_dir,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'dntmc:directory',1,dntmc_dir))
     &  dntmc_dir(1:2)='./'
      call dntmc_getspecinfo(rtdb, nmols, natom, rinit, uspc, nspc,
     &  m, SPC, tag)
C      do i=1, nmols
C        spc(i) = uspc(1) !assumes only one species
C      enddo
        if (ga_nodeid() .eq. 0) then
          do i=1, procrestart
            call dntmc_build_prcfil(prcfil, i)
            call dntmc_copy_filertdb(prcfil, rtdb, nmols, fprefix,
     &           natom, spc, dntmc_dir)
          enddo
            if(.not. rtdb_delete(rtdb, 'dntmc:procrestart'))
     &        call errquit('dntmc_restart_post:rtdb failed to delete 
     &        dntmc:procrestart'
     &        ,0,RTDB_ERR)
        endif
C End Main Program
      return
      end
C***********************************************************************
C      End of subroutine dntmc_restart_post(procrestart,rtdb)
C***********************************************************************

C***********************************************************************
      subroutine dntmc_build_prcfil(prcfil, i)
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

C***********************************************************************
C  End of subroutine dntmc_build_prcfil(prcfil, i)
C***********************************************************************

C***********************************************************************
      subroutine dntmc_copy_filertdb(prcfil, rtdb, nmols, fprefix,
     &  natom, spc, dntmc_dir)
C***********************************************************************
      Implicit none

C  Include Statements
#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"
#include "inp.fh"

C  Variable Declarations
C Everything internal
      double precision uu, usum, uavr, usqrdsum, usqrd
      double precision rsim, rsimsum, rsimavr
      double precision rmin, rmax
      double precision seed, r, numbersum, numberavr
      integer number, iacceptt, iacceptv, iacceptr, iacceptpr, iacceptpi
      integer ivtry, ittry, irtry
      integer rsmcstps, i, j, ir
      integer natom, spc, rtdb ! input
      double precision rsh
      integer nmols !input
      character*4 prcfil ! input
      character*256 fprefix ! input
      character*256 fprefixcat
      integer ncent, k, rdim
      double precision coords, dipole(3)
      dimension coords(3, MXATOM*MXMOLS)
      dimension R(MXMOLS,MXATOM,3)
      dimension NUMBER(NOB)
      dimension NUMBERSUM(NOB)
      dimension NUMBERAVR(NOB)
      dimension SEED(1)
      dimension natom(MXSP)
      dimension spc(MXMOLS)
      character*16 tag(mxatom*mxmols)
      integer ti
      character*256 dntmc_dir

C  Main Program
      write(fprefixcat,'(256(a))') (' ', i=1, 256)
      write(fprefixcat,'(5(a))')
     &      dntmc_dir(1:inp_strlen(dntmc_dir)),
     &      '/',
     &      fprefix(1:inp_strlen(fprefix)),
     &      '.restart.',
     &      prcfil(1:4)

      OPEN(UNIT = 20,FILE = fprefixcat(1:inp_strlen(fprefixcat))
     +    ,STATUS = 'OLD')

         READ(20,*)SEED(1)
         READ(20,*)UU
         READ(20,*)USUM
         READ(20,*)UAVR
         READ(20,*)USQRDSUM
         READ(20,*)USQRD
         READ(20,*)(DIPOLE(IR), IR = 1,3)
         READ(20,*)RMIN, RMAX
         READ(20,*)RSIM
         READ(20,*)RSIMSUM
         READ(20,*)RSIMAVR
         READ(20,*)IACCEPTT
         READ(20,*)IACCEPTR
         READ(20,*)IACCEPTV
         READ(20,*)IVTRY,ITTRY,IRTRY
         READ(20,*)IACCEPTPR
         READ(20,*)IACCEPTPI
         READ(20,*)rsmcstps
         ti = 1
         DO I = 1, NMOLS
            DO J = 1, NATOM(SPC(I))
               READ(20,*)tag(ti),(R(I,J,IR), IR = 1,3)
               ti = ti + 1
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
     &   call errquit('dntmc_copy_filertdb:rtdb error for seed(1)',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:uu:'//prcfil(1:4), MT_DBL, 1, UU))
     &   call errquit('dntmc_copy_filertdb:rtdb error for UU',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:usum:'//prcfil(1:4), MT_DBL, 1,
     &   usum))
     &   call errquit('dntmc_copy_filertdb:rtdb error for usum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:uavr:'//prcfil(1:4), MT_DBL, 1,
     &   uavr))
     &   call errquit('dntmc_copy_filertdb:rtdb error for uavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:usqrdsum:'//prcfil(1:4), MT_DBL, 1
     &   ,usqrdsum))
     &   call errquit('dntmc_copy_filertdb:rtdb error for usqrdsum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:usqrd:'//prcfil(1:4), MT_DBL, 1,
     &   usqrd))
     &   call errquit('dntmc_copy_filertdb:rtdb error for usqrd',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:dipole:'//prcfil(1:4), MT_DBL, 3,
     &   dipole))
     &   call errquit('dntmc_copy_filertdb:rtdb error for dipole',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rmax:'//prcfil(1:4), MT_DBL, 1,
     &   rmax))
     &   call errquit('dntmc_copy_filertdb:rtdb error for rmax',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rmin:'//prcfil(1:4), MT_DBL, 1,
     &   rmin))
     &   call errquit('dntmc_copy_filertdb:rtdb error for rmin',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsim:'//prcfil(1:4), MT_DBL, 1,
     &   rsim))
     &   call errquit('dntmc_copy_filertdb:rtdb error for rsim',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsimsum:'//prcfil(1:4), MT_DBL, 1,
     &   rsimsum))
     &   call errquit('dntmc_copy_filertdb:rtdb error for rsimsum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsimavr:'//prcfil(1:4), MT_DBL, 1,
     &   rsimavr))
     &   call errquit('dntmc_copy_filertdb:rtdb error for rsimavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptt:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptt))
     &   call errquit('dntmc_copy_filertdb:rtdb error for iacceptt',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptr:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptr))
     &   call errquit('dntmc_copy_filertdb:rtdb error for iacceptr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptv:'//prcfil(1:4), MT_INT, 1
     &   ,iacceptv))
     &   call errquit('dntmc_copy_filertdb:rtdb error for iacceptv',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:ivtry:'//prcfil(1:4), MT_INT, 1
     &   ,ivtry))
     &   call errquit('dntmc_copy_filertdb:rtdb error for ivtry',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:ittry:'//prcfil(1:4), MT_INT, 1
     &   ,ittry))
     &   call errquit('dntmc_copy_filertdb:rtdb error for ittry',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:irtry:'//prcfil(1:4), MT_INT, 1
     &   ,irtry))
     &   call errquit('dntmc_copy_filertdb:rtdb error for irtry',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptpr:'//prcfil(1:4), MT_INT,
     &   1, iacceptpr))
     &   call errquit('dntmc_copy_filertdb:rtdb error for iacceptpr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:iacceptpi:'//prcfil(1:4), MT_INT,
     &   1, iacceptpi))
     &   call errquit('dntmc_copy_filertdb:rtdb error for iacceptpi',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:number:'//prcfil(1:4), MT_INT, NOB
     &   ,number))
     &   call errquit('dntmc_copy_filertdb:rtdb error for number',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:numbersum:'//prcfil(1:4), MT_DBL,
     &   NOB,
     &   numbersum))
     &   call errquit('dntmc_copy_filertdb:rtdb error for numbersum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:numberavr:'//prcfil(1:4), MT_DBL,
     &   NOB,
     &   numberavr))
     &   call errquit('dntmc_copy_filertdb:rtdb error for numberavr',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsh:'//prcfil(1:4), MT_DBL, 1, rsh
     &   ))call errquit('dntmc_copy_filertdb:rtdb error for rsh',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:rsmcstps:'//prcfil(1:4), MT_INT, 1
     &   ,rsmcstps))
     &   call errquit('dntmc_copy_filertdb:rtdb error for rsmcstps',0,
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
     &   call errquit('dntmc_copy_filertdb:rtdb error for r',0,
     &   RTDB_ERR)

      if (.not. rtdb_cput(rtdb, 'dntmc:attag:'//prcfil(1:4),
     &   ncent, tag))
     &   call errquit('dntmc_copy_filertdb:rtdb error for attag',0,
     &   RTDB_ERR)

C  End Main Program
      return
      end

C***********************************************************************
C  End of subroutine dntmc_copy_filertdb(prcfil, rtdb, nmols, fprefix,
C     &  natom, spc, dntmc_dir)
C***********************************************************************

C***********************************************************************
      Subroutine dntmc_initialstate(nmols, rinit)
C***********************************************************************
      Implicit none
C Include Files
#include "global.fh"
#include "tcgmsg.fh"
#include "gibbsP.fh"
#include "util.fh"

C Variable Declarations
      integer nmols ! input
      double precision rinit ! input
      dimension rinit(MXMOLS,MXATOM,3)
      integer i, j, k
      double precision DIST, DELTA
      dimension delta(3)
      double precision dntmc_ran3
      logical init_condition
      double precision beta, sdbeta, cdbeta
      integer ia1, ia2, ia3
      double precision comp2, comp3, comp2r, comp3r
      double precision Maxradius(MXMOLS), radius
      integer TU, ER, IR

C Common Block Variables
      double precision r, rcm, rs ! R is output
      integer natom, spc
      integer IDUM
      double precision DISP,ALPHA
      integer ROTF
      double precision RCLCM

C Common Blocks
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)
      COMMON/SEED/IDUM
      COMMON/TRPS/DISP,ALPHA,ROTF
      COMMON/COORD2/RCLCM(3)

C Begin Subroutine
C      IDUM = NODE
C get initial structures !SPC general
      DO I = 1, nmols
        DO J=1, NATOM(SPC(I))
          R(I,J,1) = RINIT(SPC(I),J,1)
          R(I,J,2) = RINIT(SPC(I),J,2)
          R(I,J,3) = RINIT(SPC(I),J,3)
        END DO
      END DO

C Stack Molecules
      DO I = 1, nmols
        CALL dntmc_CMASS(I)
      END DO
      DO I = 1, nmols
        DO J =1, NATOM(SPC(I))
          DO K = 1,3
            R(I,J,K) = R(I,J,K) - RCM(I,K)
          END DO
        END DO
      END DO

C Find Maximum Extent +10%
C Center Molecule
      CALL dntmc_CLUSTERCMASS(NMOLS)

      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO
C Find Maximum Extent +10%
      DO I =1, NMOLS
        CALL dntmc_CMASS(I)
      END DO
      DO I = 1, NMOLS
        Maxradius(I) = 0.D+00
        DO J = 1, NATOM(SPC(I))
          radius = 0.00D+00
          DO K = 1, 3
            radius = radius + (R(I,J,K) - RCM(I,K))**2.D+00
          END DO
          radius = DSQRT(radius)
          if (radius .gt. Maxradius(I)) Maxradius(I) = radius
        END DO
        Maxradius(I) = Maxradius(I) * 1.10D+00
        if(util_print('MCout',print_debug)) then
        if (ga_nodeid() .eq. 0) 
     &  WRITE(12,*)'Molecule ',I,' has an ',Maxradius(I),
     &  ' Angstrom radial extent.'
        endif
      END DO

      DO TU = 1, NMOLS ! Selects molecules incrementally
      init_condition = .false.
 1357 CONTINUE ! This is my while loop
C Translation
      I = TU
C      DO I = 1, NMOLS
C     PICK DISPLACEMENTS ALONG X,Y,Z ! 10 times larger step size
         DELTA(1) = DISP*10.D+00*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
         DELTA(2) = DISP*10.D+00*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
         DELTA(3) = DISP*10.D+00*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) + DELTA(1)
            R(I,J,2) = R(I,J,2) + DELTA(2)
            R(I,J,3) = R(I,J,3) + DELTA(3)
         END DO
C      END DO

C Center Molecule
      CALL dntmc_CLUSTERCMASS(NMOLS)

      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

C Rotation
      I = TU
C      DO I = 1, NMOLS
C     PICK AN ANGLE ! 10 times larger rotational angle
      BETA = ALPHA*10.D+00*(2.D+00*dntmc_RAN3(IDUM) - 1.D+00)
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
C      END DO

C Center Molecule
      CALL dntmc_CLUSTERCMASS(NMOLS)

      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

C Check for uniformity pairwise
      ER = 0
      DO I = 1, NMOLS
        CALL dntmc_CMASS(I)
      END DO
      I = TU
C      DO I = 1, NMOLS - 1
        DO J = 1, NMOLS
          if (J .ne. I) then
            radius = 0.D+00
            DO K = 1, 3
              radius = radius + (RCM(I,K) - RCM(J,K))**2.D+00
            END DO
            radius = DSQRT(radius)
            if ((Maxradius(I) + Maxradius(J)) .gt. radius) ER = ER+1
          endif
        END DO
C      END DO
      if (ER .eq. 0) init_condition = .true.
      if (.not. init_condition) goto 1357 ! End of While Loop
      if(util_print('MCout',print_debug)) then
      if (ga_nodeid() .eq. 0) 
     &    WRITE(12,*)'Molecule ',I, ' initial position found.'
      endif
      
      ENDDO ! TU Loop

      if (ga_nodeid() .eq. 0) then
        if(util_print('MCout',print_debug)) then
        Write(12,*) 'Initial Coordinates after Generation are:'
        DO I = 1, NMOLS
          DO J = 1, NATOM(SPC(I))
            WRITE(12,*)(R(I,J,IR),IR=1,3)
          ENDDO
        ENDDO
        CALL FLUSH(12)
        endif
      endif

C End Subroutine
      return
      end

C***********************************************************************
C  End of Subroutine dntmc_initialstate(nmols, rinit)
C***********************************************************************
C***********************************************************************
      Subroutine dntmc_initialstate2(nmols, rinit, radiuscube)
C***********************************************************************
      Implicit none
C Include Files
#include "global.fh"
#include "tcgmsg.fh"
#include "gibbsP.fh"
#include "util.fh"

C Variable Declarations
      double precision radiuscube ! input
      integer nmols ! input
      double precision rinit !input
      dimension rinit(MXMOLS,MXATOM,3)
      integer i, j, k, P
      double precision DELTA
      dimension delta(3)
      double precision dntmc_ran3
      logical overlapflag
      double precision beta, sdbeta, cdbeta
      integer ia1, ia2, ia3
      double precision comp2, comp3, comp2r, comp3r
      double precision Maxradius(MXMOLS), radius
      integer IR
      double precision dr
 

C Common Block Variables
      double precision r, rcm, rs ! R is output
      integer natom, spc
      integer IDUM
      integer ROTF
      double precision RCLCM
      double precision e,ki,esqrd,kb,pi

C Common Blocks
      COMMON/COORD/R(MXMOLS,MXATOM,3),RCM(MXMOLS,3),
     +RS(MXMOLS,MXATOM,3)
      COMMON/NNATOMS/NATOM(MXSP)
      COMMON/SPECIES/SPC(MXMOLS)
      COMMON/SEED/IDUM
      COMMON/COORD2/RCLCM(3)
      COMMON/CONSTANTS/e,ki,esqrd,kb,PI

C Begin Subroutine
C Begin Placement loop
      Do I = 1, NMOLS

C Place monomer
      Do J = 1, NATOM(SPC(I))
        Do K = 1,3
          R(I,J,K) = RINIT(SPC(I),J,K)
        EndDo
      EndDo
C Center at orgin
      CALL dntmc_CMASS(I)
      Do J = 1, NATOM(SPC(I))
        Do K = 1,3
          R(I,J,K) = R(I,J,K) - RCM(I,K)
        EndDo
      ENDDO

C Define non-overlapping sphere + 10%
      CALL dntmc_CMASS(I)
      Maxradius(I) = 0.00D+00
      Do J = 1, NATOM(SPC(I))
        radius = 0.00D+00
        Do K = 1,3
          radius = radius + (R(I,J,K) - RCM(I,K))**2.0D+00
        ENDDO
        radius = DSQRT(radius)
        if (radius .gt. Maxradius(I)) Maxradius(I) = radius
      ENDDO
      Maxradius(I) = Maxradius(I) * 1.10D+00
        if(util_print('MCout',print_debug)) then
        if (ga_nodeid() .eq. 0)
     &  WRITE(12,*)'Molecule ',I,' has an ',Maxradius(I),
     &  ' Angstrom radial extent.'
        endif

C Check to see if radius cube is large enough
      radius = 0.00D+00
      Do J = 1, I
        radius = radius + (8.00D+00*(Maxradius(J)**3.00D+00))
      EndDo
      if (radius .gt. (8.00D+00/3.00D+00/DSQRT(3.00D+00)*
     &   (radiuscube**3.00D+00)))
     & then
        radiuscube = 
     &  (radius*3.00D+00*DSQRT(3.00D+00)/8.00D+00)**(1.00D+00/3.00D+00)
      endif

      if (I .eq. 1) goto 5678
C Radomly place in cube of given radius
      P = 0
 1234 Continue
      Do K = 1,3
        Delta(K) = radiuscube/DSQRT(3.00D+00)*
     &   ((2.0D+00*dntmc_RAN3(IDUM)) - 1.0D+00)
      ENDDO
      radius = 0.00D+00
      Do K = 1,3
        radius = radius + (Delta(K)**2.00D+00)
      ENDDO
      radius = DSQRT(radius)
      if (radius .gt. radiuscube) goto 1234
      Do J = 1, I-1
        Call dntmc_CMASS(J)
      ENDDO
      overlapflag = .false.
      Do J = 1, I-1
        radius = 0.00D+00
        Do K = 1,3
          radius = radius +((RCM(J,K) - Delta(K))**2.00D+00)
        EndDo
        radius = DSQRT(radius)
        if (radius .lt. (Maxradius(I) + Maxradius(J))) then
          overlapflag = .true.
        endif
      EndDo
C    if having a little trouble increase the radius 10%
      if (overlapflag) then
        if (P .eq. NOB) then
          radiuscube = radiuscube * 1.10D+00
          P = 0
        endif
        P = P + 1
        goto 1234
      endif
C Placement worked now moving monomer
      Do J = 1, NATOM(SPC(I))
        Do K = 1,3
          R(I,J,K) = R(I,J,K) + DELTA(K)
        EndDo
      EndDo
C Check that no monomer has fallen out of volume
      Call dntmc_ClusterCMASS(I)
      Do J = 1, I
          DR = 0.00D+00
        Do K = 1, 3
          DR = DR + (RCM(J,K)-RCLCM(K))**2.00D+00
        ENDDO
        if (DR .gt. radiuscube**2.00D+00) then
          Do P = 1, NATOM(SPC(I))
          Do K = 1,3
            R(I,P,K) = R(I,P,K) - DELTA(K)
          EndDo
          ENDDO
          P = 1
          goto 1234
        endif
      EndDo

C Recenter cluster to this point at COM.
      Call dntmc_ClusterCMASS(I)
      Do J = 1, I
        Do K = 1, NATOM(SPC(J))
          Do P = 1,3
            R(J,K,P) = R(J,K,P) - RCLCM(P)
          EndDo
        EndDo
      EndDo
        
 5678 CONTINUE
      End DO ! I Loop

C Rotate Monomers along all directions
      Do I = 1, NMOLS
        Do K = 1, 3
        Beta = 2.00D+00*PI*((2.0D+00*dntmc_Ran3(IDUM)) - 1.00D+00)
        SDBETA = DSIN(BETA)
        CDBETA = DCOS(BETA)
        IA1 = K
        IA2 = IA1 + 1
        IA3 = IA1 + 2
        if(IA2 .gt. 3) IA2 = IA2 - 3
        if(IA3 .gt. 3) IA3 = IA3 - 3
        CAll dntmc_cmass(I)
          Do J = 1, NATOM(SPC(I))
            COMP2 = (R(I,J,IA2) - RCM(I,IA2))
            COMP3 = (R(I,J,IA3) - RCM(I,IA3))
            COMP2R = CDBETA*COMP2+SDBETA*COMP3
            COMP3R = -SDBETA*COMP2+CDBETA*COMP3
            R(I,J,IA2) = RCM(I,IA2) + COMP2R
            R(I,J,IA3) = RCM(I,IA3) + COMP3R
          EndDo
        EndDo
      EndDo

C Placement of monomers is complete
      if (ga_nodeid() .eq. 0) then
        if(util_print('MCout',print_debug)) then
        Write(12,*) 'Initial Coordinates after Generation are:'
        DO I = 1, NMOLS
          DO J = 1, NATOM(SPC(I))
            WRITE(12,*)(R(I,J,IR),IR=1,3)
          ENDDO
        ENDDO
        CALL FLUSH(12)
        endif
      endif

C End of subroutine
      return
      end

C***********************************************************************
C  End of Subroutine dntmc_initialstate2(nmols, rinit, radiuscube)
C***********************************************************************

