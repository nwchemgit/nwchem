C$Id$
C23456789012345678901234567890123456789012345678901234567890123456789012
C     modified on oct 3, 1999 for restart on mpp1
C     modified on oct 25, 1999 for random number generator restart
C     modified on may 10, 2004 for colony2a
C     modified on april 15, 2005 for ran3
C     modified on and about april 2, 2007 to clean up code (lcrosby)
C     version 5.0  Dec. 20, 2007 full release version (lcrosby)
C     modified on June 6, 2008 to include performace timers (lcrosby)
C     modified on June 10, 2008 to reorganize code (lcrosby)
C     modified on June 11, 2008 to support multiple species (lcrosby)


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
#include "pstat.fh" 

C#include "errquit.fh" ! for NWChem errquit function
C#include "rtdb.fh" ! for rtdb operations
C#include "mafdecls.fh" ! for MA data types
C#include "tcgmsg.fh" ! for tcgmsg support
C#include "global.fh" ! for GA support
C#include "msgids.fh" ! for message IDs
C#include "msgtypesf.h" ! for data types for message IDs
C#include "gibbsP.fh" ! for gibbsP header file
C#include "util.fh" ! for NWChem utilites
C#include "inp.fh" ! for input processing
C#include "util_sgroup.fh" ! for subgroup support
C#include "pstat.fh" ! for performance timers



C Variable declarations
      integer rtdb
      integer myid, idbig, inodesbig
      integer idmedium, idzero
      integer g_pr
      integer iacceptt, iacceptr, iacceptv, iacceptpr, iacceptpi
      integer ivtry, ittry, irtry
      integer cont1, cont2, mprnt
      double precision temp, uu, uavr, usqrd, rsimsum, rsimavr
      logical rconf
      integer lmin, lmax
      character*4 prcfil
      character*4 prcfil2
      integer procrestart
      character*256 fprefix
      double precision rest
      logical flag
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

C Array Variables
      integer l_groupnums, k_groupnums
      integer l_nda_old, k_nda_old
      integer l_nda, k_nda
      integer l_ndsig, k_ndsig
      integer l_testme, k_testme
      integer l_rb, k_rb
      integer nob, numgroups
    
      integer l_natom, k_natom
      integer l_rinit, k_rinit
      integer l_uspc, k_uspc
      integer l_m, k_m
      integer l_spc, k_spc
      integer l_tag, k_tag
      integer nmols, nspc, nattotal, natmax

      integer l_mtot, k_mtot
      integer l_r, k_r
      integer l_rcm, k_rcm
      integer l_rs, k_rs
      integer l_ri, k_ri

      integer l_number, k_number
      integer l_numbersum, k_numbersum 
      integer l_numberavr, k_numberavr

      integer l_spctag, k_spctag

      double precision dntmc_get_array_element
      double precision dntmc_get_array_element2

C Variable declarations for common block
      integer L

C Timer Common block variables
      integer i_mcloop, i_together
      character*16 junk
      logical statuspstat

C Common Blocks
      COMMON/MSTP/L

C Common handles
      COMMON/pstathandels/i_mcloop, i_together, junk

C Main program
      dntmc_gibbs = .false. ! function value

C get print levels from input
      call util_print_push
      call util_print_rtdb_load(rtdb,'dntmc')

C allocate arrays
      call dntmc_allocate_arrays(l_groupnums, k_groupnums,
     & l_nda_old, k_nda_old, l_nda, k_nda, l_ndsig, k_ndsig, l_testme,
     & k_testme, l_rb, k_rb, nob, numgroups, l_natom, k_natom, l_rinit,
     & k_rinit, l_uspc, k_uspc, l_m, k_m, l_spc, k_spc, l_tag, k_tag,
     & nmols, nspc, nattotal, natmax, l_mtot, k_mtot, l_r, k_r, l_rcm,
     & k_rcm, l_rs, k_rs, l_ri, k_ri, rtdb, l_number, k_number,
     & l_numbersum, k_numbersum, l_numberavr, k_numberavr, l_spctag,
     & k_spctag)

C initial print
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
      if (.not. rtdb_get(rtdb, 'dntmc:nob', MT_INT, 1, nob))
     &  call errquit
     &('gibbs: failed to get nob from rtdb',0,RTDB_ERR)
        call dntmc_restart_post(procrestart,rtdb,nob)
      goto 9898 ! Speed to end of code
      endif

C Get Info on processes
      myid = ga_nodeid()
      idbig = ga_pgroup_get_default()
      inodesbig = ga_nnodes()

C Setup Timers
      if(util_print('Timers',print_debug)) then
      statuspstat = Pstat_Init(2,0,junk)
      statuspstat =
     &  pstat_allocate('dntmc:mcloop',Pstat_QStat,0,junk,i_mcloop)
      endif

C Setup GA
      if(ga_nodeid() .eq. 0) then
        WRITE(6,*)'******** Splitting Into SubGroups *********'
        WRITE(6,*)'*******************************************'
        call flush(6)
      endif

      call dntmc_setup_ga
     &(rtdb, myid, idbig, inodesbig,
     &idmedium, idzero, g_pr, int_mb(k_groupnums),numgroups,nob,nspc) ! Now in Subgroups

C Initialize group related information
      if(ga_nodeid() .eq. 0) then
      call dntmc_setup_group
     &(rtdb, idmedium, myid, int_mb(k_groupnums),prcfil, prcfil2,
     & fprefix, numgroups)
      else
      call dntmc_write_prcfil(idmedium, int_mb(k_groupnums), prcfil,
     &   numgroups)
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

      call dntmc_setup_dntmc(rtdb, iacceptt, iacceptr, iacceptv, ivtry,
     &  ittry, irtry, iacceptpr, iacceptpi, cont1, cont2, nmols, mprnt,
     &  temp, rconf, lmin, lmax, myid, idmedium, int_mb(k_groupnums),
     &  prcfil, uu, uavr, usqrd, rsimsum, rsimavr, nspc, dipole,
     &  conv_perc, numgroups, nattotal, natmax, int_mb(k_natom),
     &  dbl_mb(k_rinit), int_mb(k_uspc), dbl_mb(k_m), int_mb(k_spc),
     &  byte_mb(k_tag), dbl_mb(k_mtot), dbl_mb(k_r), dbl_mb(k_rcm),
     &  dbl_mb(k_rs), dbl_mb(k_ri), nob, int_mb(k_number),
     *  dbl_mb(k_numbersum), dbl_mb(k_numberavr), int_mb(k_spctag))

C Starting Monte Carlo Loop
      if(ga_pgroup_nodeid(idbig) .eq. 0) then
        WRITE(6,*)'********** Starting MC Loop ***************'
        WRITE(6,*)'*******************************************'
        call flush(6)
      endif

      loop_avg = 0.00D+00
      DO j = 1, nspc
       DO i = 1, NOB
        call dntmc_put_array_element2(dbl_mb(k_nda_old),nob,nspc,
     &     i,j,0.00D+00)
       ENDDO
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
        endif
      endif
C     We are starting loops, so flush all output buffers
      call flush(6)

C start Timers
      if(util_print('Timers',print_debug)) then
      Call pstat_on(i_mcloop)
      endif

C MC Do Loop exits upon convergence (flag) or L > LMAX
      Do L = LMIN, LMAX
      loop_start = util_wallsec() ! for time checking

      call dntmc_mcloop(rtdb,temp,UU,nmols,nspc,iacceptt,iacceptr
     &  ,iacceptv,iacceptpr,iacceptpi,uavr,rsimavr,rconf,prcfil,fprefix
     &  ,mprnt,usqrd,rsimsum,idmedium,idzero,cont1,g_pr,rest,flag
     &  ,dbl_mb(k_nda_old),conv_perc,ivtry,ittry,irtry,dipole,idbig,nob,
     & numgroups,natmax,dbl_mb(k_rb),dbl_mb(k_testme),
     & dbl_mb(k_numbersum),dbl_mb(k_numberavr),dbl_mb(k_nda),
     & dbl_mb(k_ndsig),dbl_mb(k_m),dbl_mb(k_mtot),int_mb(k_uspc),
     & nattotal,dbl_mb(k_r),dbl_mb(k_rcm),dbl_mb(k_rs),int_mb(k_natom),
     & int_mb(k_spc),byte_mb(k_tag),int_mb(k_number),int_mb(k_spctag))

C     Time Checking
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
      enddo ! End of DO L Loop
      lnow = l

C Stop Timers
      if(util_print('Timers',print_debug)) then
      Call pstat_off(i_mcloop)
      endif
C Act like not converged if run out of iterations 
      goto 299 ! Not Converged

C End of Monte Carlo Loop 1st.
 199  CONTINUE ! Converged
C Finish up MC Loop
      IF(CONT1.EQ.1)THEN
      if (ga_nodeid() .eq. 0) then
        call dntmc_writefiles(uu, uavr, usqrd, ivtry, ittry, irtry,
     +    rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +    iacceptpi, rsh, prcfil, nmols, lnow, fprefix, dipole, nob,
     &  natmax,nspc, nattotal, dbl_mb(k_r), dbl_mb(k_rcm), dbl_mb(k_rs),
     &  int_mb(k_number), dbl_mb(k_numbersum), dbl_mb(k_numberavr),
     &  int_mb(k_natom), int_mb(k_spc), byte_mb(k_tag),int_mb(k_spctag))
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
      WRITE(3,*)((L-1)*util_sgroup_numgroups()), 'Total Configurations'
      DO J=1, nspc
      WRITE(3,*)'Species number ',J
      DO I=1, NOB
      WRITE(3,*)dntmc_get_array_element(dbl_mb(k_rb),nob,I),
     & dntmc_get_array_element2(dbl_mb(k_nda),nob,nspc,I,J),
     & dntmc_get_array_element2(dbl_mb(k_ndsig),nob,nspc,i,j)
      ENDDO
      ENDDO
      WRITE(3,*)'***************************************'
      call flush(3)
      endif
      endif

C Job Exiting 
      goto 999 ! End Job

C End of Monte Carlo Loop 2nd
 299  CONTINUE ! Not Converged

C Finish up Job (Not Converged)
      IF(CONT1.EQ.1)THEN
      if (ga_nodeid() .eq. 0) then
        call dntmc_writefiles(uu, uavr, usqrd, ivtry, ittry, irtry,
     +    rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +    iacceptpi, rsh, prcfil, nmols, lnow, fprefix, dipole, nob,
     &  natmax,nspc, nattotal, dbl_mb(k_r), dbl_mb(k_rcm), dbl_mb(k_rs),
     &  int_mb(k_number), dbl_mb(k_numbersum), dbl_mb(k_numberavr),
     &  int_mb(k_natom), int_mb(k_spc), byte_mb(k_tag),int_mb(k_spctag))
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
      DO j=1, nspc
      Write(3,*)'Species number ',J
      DO I=1, NOB
      WRITE(3,*)dntmc_get_array_element(dbl_mb(k_rb),nob,I),
     & dntmc_get_array_element2(dbl_mb(k_nda),nob,nspc,I,J),
     & dntmc_get_array_element2(dbl_mb(k_ndsig),nob,nspc,i,j)
      ENDDO
      ENDDO
      WRITE(3,*)'***************************************'
      call flush(3)
      endif
      endif

C End of DNTMC Run
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

C End Performace Timers
      if(util_print('Timers',print_debug)) then
        statuspstat = pstat_terminate()
      endif

C Clean up arrays
      call dntmc_deallocate_arrays(l_groupnums)

C Print ending messages
      call util_print_pop

C Very End of RUN
 9898 dntmc_gibbs = .true.
      return
      end
C***********************************************************************
C***********************************************************************
      Subroutine dntmc_allocate_arrays(l_groupnums, k_groupnums,
     & l_nda_old, k_nda_old, l_nda, k_nda, l_ndsig, k_ndsig, l_testme,
     & k_testme, l_rb, k_rb, nob, numgroups, l_natom, k_natom, l_rinit,
     & k_rinit, l_uspc, k_uspc, l_m, k_m, l_spc, k_spc, l_tag, k_tag,
     & nmols, nspc, nattotal, natmax, l_mtot, k_mtot, l_r, k_r, l_rcm,
     & k_rcm, l_rs, k_rs, l_ri, k_ri, rtdb, l_number, k_number,
     & l_numbersum, k_numbersum, l_numberavr, k_numberavr, l_spctag,
     & k_spctag)
C***********************************************************************
      Implicit none

C Include Statements
#include "rtdb.fh" 
#include "mafdecls.fh" 
#include "errquit.fh" 

C#include "rtdb.fh" !RTDB routines
C#include "mafdecls.fh" !MA allocator
C#include "errquit.fh" !error quit

C Variables
      integer rtdb

      integer l_groupnums, k_groupnums
      integer l_nda_old, k_nda_old
      integer l_nda, k_nda
      integer l_ndsig, k_ndsig
      integer l_testme, k_testme
      integer l_rb, k_rb
      integer nob, numgroups

      integer l_natom, k_natom
      integer l_rinit, k_rinit
      integer l_uspc, k_uspc
      integer l_m, k_m
      integer l_spc, k_spc
      integer l_tag, k_tag
      integer nmols, nspc, nattotal, natmax

      integer l_mtot, k_mtot
      integer l_r, k_r
      integer l_rcm, k_rcm
      integer l_rs, k_rs
      integer l_ri, k_ri

      integer l_number, k_number
      integer l_numbersum, k_numbersum
      integer l_numberavr, k_numberavr

      integer l_spctag, k_spctag

C Subroutine
C Get some lengths
      if (.not. rtdb_get(rtdb, 'dntmc:nob', MT_INT, 1, nob))
     &  call errquit
     &('dntmc_allocate_arrays: failed to get nob from rtdb',0,RTDB_ERR)

      if (.not. rtdb_get(rtdb,'subgroups_number', MT_INT, 1, numgroups))
     & numgroups = 1

C Allocate memory
      if (.not. ma_push_get(MT_INT, numgroups, 'dntmc: groupnums',
     & l_groupnums, k_groupnums)) call errquit
     &('dntmc_allocate_arrays: cannot allocate groupnums',0,MA_ERR)

      call dntmc_getspecinfo(rtdb, nmols, l_natom, k_natom, l_rinit,
     & k_rinit, l_uspc, k_uspc, nspc, l_m, k_m, l_SPC, k_spc, l_tag
     &, k_tag, nattotal, natmax)

C      call test_arrays(nmols,nspc,nattotal,natmax,int_mb(k_natom)
C     & ,dbl_mb(k_rinit),int_mb(k_uspc),dbl_mb(k_m),int_mb(k_spc),
C     & byte_mb(k_tag))

C Changed to dimentions (nob,nspc)
      if (.not. ma_push_get(MT_DBL, nob*nspc, 'dntmc: nda_old',
     & l_nda_old, k_nda_old)) call errquit
     &('dntmc_allocate_arrays: cannot allocate nda_old',0,MA_ERR)

      if (.not. ma_push_get(MT_DBL, nob*nspc, 'dntmc: nda',
     & l_nda, k_nda)) call errquit
     &('dntmc_allocate_arrays: cannot allocate nda',0,MA_ERR)

      if (.not. ma_push_get(MT_DBL, nob*nspc, 'dntmc: ndsig',
     & l_ndsig, k_ndsig)) call errquit
     &('dntmc_allocate_arrays: cannot allocate ndsig',0,MA_ERR)

C Changed to dimentions (nob,nspc,numgroups)
      if (.not. ma_push_get(MT_DBL, nob*numgroups*nspc, 'dntmc: testme',
     & l_testme, k_testme)) call errquit
     &('dntmc_allocate_arrays: cannot allocate testme',0,MA_ERR)

C End of Changes

      if (.not. ma_push_get(MT_DBL, nob, 'dntmc: rb',
     & l_rb, k_rb)) call errquit
     &('dntmc_allocate_arrays: cannot allocate rb',0,MA_ERR)

      if (.not. ma_push_get(MT_DBL, nspc, 'dntmc: mtot',
     & l_mtot, k_mtot)) call errquit
     &('dntmc_allocate_arrays: cannot allocate mtot',0,MA_ERR)

      if (.not. ma_push_get(MT_DBL, nmols*natmax*3, 'dntmc: r',
     & l_r, k_r)) call errquit
     &('dntmc_allocate_arrays: cannot allocate r',0,MA_ERR)

      if (.not. ma_push_get(MT_DBL, nmols*3, 'dntmc: rcm',
     & l_rcm, k_rcm)) call errquit
     &('dntmc_allocate_arrays: cannot allocate rcm',0,MA_ERR)

      if (.not. ma_push_get(MT_DBL, nmols*natmax*3, 'dntmc: rs',
     & l_rs, k_rs)) call errquit
     &('dntmc_allocate_arrays: cannot allocate rs',0,MA_ERR)

      if (.not. ma_push_get(MT_DBL, nmols*natmax*3, 'dntmc: ri',
     & l_ri, k_ri)) call errquit
     &('dntmc_allocate_arrays: cannot allocate ri',0,MA_ERR)

C Changed to dimentions (nob, nspc)
      if (.not. ma_push_get(MT_INT, nob*nspc, 'dntmc: number',
     & l_number, k_number)) call errquit
     &('dntmc_allocate_arrays: cannot allocate number',0,MA_ERR)

      if (.not. ma_push_get(MT_DBL, nob*nspc, 'dntmc: numbersum',
     & l_numbersum, k_numbersum)) call errquit
     &('dntmc_allocate_arrays: cannot allocate numbersum',0,MA_ERR)

      if (.not. ma_push_get(MT_DBL, nob*nspc, 'dntmc: numberavr',
     & l_numberavr, k_numberavr)) call errquit
     &('dntmc_allocate_arrays: cannot allocate numberavr',0,MA_ERR)

      if (.not. ma_push_get(MT_INT, nspc, 'dntmc: spctag',
     & l_spctag, k_spctag)) call errquit
     &('dntmc_allocate_arrays: cannot allocate spctag',0,MA_ERR)


C End of Changes

C End of Routine
      return
      end
C***********************************************************************
C End of Subroutine dntmc_allocate_arrays(l_groupnums, k_groupnums,
C     & l_nda_old, k_nda_old, l_nda, k_nda, l_ndsig, k_ndsig, l_testme,
C     & k_testme, l_rb, k_rb, nob, numgroups, l_natom, k_natom, l_rinit,
C     & k_rinit, l_uspc, k_uspc, l_m, k_m, l_spc, k_spc, l_tag, k_tag,
C     & nmols, nspc, nattotal, natmax, l_mtot, k_mtot, l_r, k_r, l_rcm,
C     & k_rcm, l_rs, k_rs, l_ri, k_ri, rtdb, l_number, k_number,
C     & l_numbersum, k_numbersum, l_numberavr, k_numberavr, l_spctag,
C     & k_spctag)
C***********************************************************************
C***********************************************************************
      Subroutine dntmc_deallocate_arrays(l_groupnums)
C***********************************************************************
      Implicit none

C Include Statements
#include "mafdecls.fh" 
#include "errquit.fh" 

C#include "mafdecls.fh" !MA allocator
C#include "errquit.fh" !error quit

C Variables
      integer l_groupnums

C Free memory by stack chop
      if (.not. ma_chop_stack(l_groupnums))
     & call errquit('dntmc_deallocate_arrays:chop of l_groupnums failed'
     &,0,MA_ERR)

C End of routine
      return
      end
C***********************************************************************
C End of Subroutine dntmc_deallocate_arrays(l_groupnums)
C***********************************************************************
C***********************************************************************
      double precision function 
     &      dntmc_get_array_element(k_tag,length,element)
C***********************************************************************
      Implicit none

C Variable Declaration
      integer length, element
      double precision k_tag(length)

C Function
      dntmc_get_array_element = k_tag(element)

C End of Function
      return
      end
C***********************************************************************
C End of double precision function
C     &      dntmc_get_array_element(k_tag,length,element)
C***********************************************************************
C***********************************************************************
      double precision function
     &      dntmc_get_array_element2(k_tag,length,length2,element,
     &      element2)
C***********************************************************************
      Implicit none

C Variable Declaration
      integer length, element
      integer length2, element2
      double precision k_tag(length,length2)

C Function
      dntmc_get_array_element2 = k_tag(element,element2)

C End of Function
      return
      end
C***********************************************************************
C End of double precision function
C     &      dntmc_get_array_element2(k_tag,length,element)
C***********************************************************************
C***********************************************************************
      integer function
     &      dntmc_get_iarray_element(k_tag,length,element)
C***********************************************************************
      Implicit none

C Variable Declaration
      integer length, element
      integer k_tag(length)

C Function
      dntmc_get_iarray_element = k_tag(element)

C End of Function
      return
      end
C***********************************************************************
C End of integer function
C     &      dntmc_get_iarray_element(k_tag,length,element)
C***********************************************************************
C***********************************************************************
      Subroutine dntmc_put_array_element(k_tag,length,element,data1)
C***********************************************************************
      Implicit none

C Variable Declaration
      integer length, element
      double precision k_tag(length)
      double precision data1

C Routine
      k_tag(element) = data1

C End of Routine
      return
      end
C***********************************************************************
C End of Subroutine dntmc_put_array_element(k_tag,length,element,data1)
C***********************************************************************
C***********************************************************************
      Subroutine dntmc_put_array_element2(k_tag,length,length2,element
     &  ,element2,data1)
C***********************************************************************
      Implicit none

C Variable Declaration
      integer length, element
      integer length2, element2
      double precision k_tag(length,length2)
      double precision data1

C Routine
      k_tag(element,element2) = data1

C End of Routine
      return
      end
C***********************************************************************
C End of Subroutine dntmc_put_array_element2(k_tag,length,length2,element
C     ,element2,data1)
C***********************************************************************
C***********************************************************************
      Subroutine dntmc_put_iarray_element(k_tag,length,element,data1)
C***********************************************************************
      Implicit none

C Variable Declaration
      integer length, element
      integer k_tag(length)
      integer data1

C Routine
      k_tag(element) = data1

C End of Routine
      return
      end
C***********************************************************************
C End of Subroutine dntmc_put_iarray_element(k_tag,length,element,data1)
C***********************************************************************
C***********************************************************************
      subroutine dntmc_restart_post(procrestart,rtdb,nob)
C***********************************************************************
      Implicit None

C Include Statements
#include "tcgmsg.fh"
#include "global.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"

C Variable Declarations
      integer rtdb, procrestart ! input
      character*256 fprefix
      integer nmols, nspc
      integer i 
      character*4 prcfil
      character*256 dntmc_dir
      integer nattotal, natmax, nob

      integer l_rinit, k_rinit
      integer l_uspc, k_uspc
      integer l_m, k_m
      integer l_spc, k_spc
      integer l_tag, k_tag
      integer l_natom, k_natom

C Everything else is internal

C Main Program
      write(fprefix,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'file_prefix',1,fprefix))
     &  call errquit('dntmc_restart_post:rtdb get file_prefix failed',0
     &  ,RTDB_ERR)
      write(dntmc_dir,'(256(a))') (' ', i=1,256)
      if(.not. rtdb_cget(rtdb,'dntmc:directory',1,dntmc_dir))
     &  dntmc_dir(1:2)='./'

      call dntmc_getspecinfo(rtdb, nmols, l_natom, k_natom, l_rinit,
     & k_rinit, l_uspc, k_uspc, nspc, l_m, k_m, l_SPC, k_spc, l_tag
     &, k_tag, nattotal, natmax)

        if (ga_nodeid() .eq. 0) then
          do i=1, procrestart
            call dntmc_build_prcfil(prcfil, i)
            call dntmc_copy_filertdb(prcfil, rtdb, nmols, fprefix,
     &           int_mb(k_natom), int_mb(k_spc), dntmc_dir,
     &           nspc, nattotal, natmax, nob)
          enddo
            if(.not. rtdb_delete(rtdb, 'dntmc:procrestart'))
     &        call errquit('dntmc_restart_post:rtdb failed to delete
     &        dntmc:procrestart'
     &        ,0,RTDB_ERR)
        endif

C Sync processes
        call ga_sync()

C Clean up memory from dntmc_getspecinfo
      if (.not.ma_chop_stack(l_natom)) call errquit
     &    ('dntmc_restart_post: ma chop failed for natom', 1, MA_ERR)

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
     &  natom, spc, dntmc_dir, nspc, nattotal, natmax, nob)
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
      integer nspc, nattotal, natmax
      double precision uu, usum, uavr, usqrdsum, usqrd
      double precision rsim, rsimsum, rsimavr
      double precision rmin, rmax
      double precision  r, numbersum, numberavr
      integer number, iacceptt, iacceptv, iacceptr, iacceptpr, iacceptpi
      integer ivtry, ittry, irtry
      integer rsmcstps, i, j, ir
      integer rtdb ! input
      double precision rsh
      integer nmols !input
      character*4 prcfil ! input
      character*256 fprefix ! input
      character*256 fprefixcat
      integer ncent, k, rdim
      double precision coords, dipole(3)
      double precision dntSEED(1)
      integer ti
      character*256 dntmc_dir

      integer nob
      integer natom(nspc), spc(nmols)
      dimension coords(3, nattotal)
      dimension R(nmols,natmax,3)

      integer spctag(nspc)

C Changed Dimensions to (nob, nspc)
      dimension NUMBER(NOB,nspc)
      dimension NUMBERSUM(NOB,nspc)
      dimension NUMBERAVR(NOB,nspc)
      character*16 tag(nattotal)

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

         READ(20,*)dntSEED(1)
         READ(20,*)UU
         READ(20,*)USUM
         READ(20,*)UAVR
         READ(20,*)USQRDSUM
         READ(20,*)USQRD
         READ(20,*)(DIPOLE(IR), IR = 1,3)
         READ(20,*)RMIN, RMAX
         READ(20,*)RSIM
         READ(20,*)(spctag(I),I=1,nspc)
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

C Changed format to columns for each spc
         DO I = 1, NOB
            READ(20,*)(NUMBER(I,J),J=1,nspc)
            READ(20,*)(NUMBERSUM(I,J),J=1,nspc)
            READ(20,*)RSH,(NUMBERAVR(I,J),J=1,nspc)
         ENDDO

         CLOSE(20)

      if (.not. rtdb_put(rtdb, 'dntmc:seed:'//prcfil(1:4), MT_DBL, 1,
     &   dntSEED(1)))
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

C Change dimenstions to (nob,nspc)
      if (.not. rtdb_put(rtdb, 'dntmc:number:'//prcfil(1:4), MT_INT,
     &   NOB*nspc,number))
     &   call errquit('dntmc_copy_filertdb:rtdb error for number',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:numbersum:'//prcfil(1:4), MT_DBL,
     &   NOB*nspc,
     &   numbersum))
     &   call errquit('dntmc_copy_filertdb:rtdb error for numbersum',0,
     &   RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'dntmc:numberavr:'//prcfil(1:4), MT_DBL,
     &   NOB*nspc,
     &   numberavr))
     &   call errquit('dntmc_copy_filertdb:rtdb error for numberavr',0,
     &   RTDB_ERR)
C End Change

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

      if (.not. rtdb_put(rtdb, 'dntmc:spctag:'//prcfil(1:4), MT_INT,
     &   nspc,spctag))
     &   call errquit('dntmc_copy_filertdb:rtdb error for spctag',0,
     &   RTDB_ERR)


C  End Main Program
      return
      end

C***********************************************************************
C  End of subroutine dntmc_copy_filertdb(prcfil, rtdb, nmols, fprefix,
C     &  natom, spc, dntmc_dir)
C***********************************************************************
C***********************************************************************
      subroutine dntmc_setup_ga
     &(rtdb, myid, idbig, inodesbig,
     &idmedium, idzero, g_pr, groupnums, numgroups, nob, nspc)
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
C timers
#include "pstat.fh"

C Variable Declarations
      integer rtdb ! input
      integer nob, nspc
      integer myid, idbig, inodesbig, numgroups ! input
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
      dimension procnums(numgroups)
      dimension groupnums(numgroups)
      dimension bigproc(inodesbig,1)

C timers common
      integer i_mcloop,i_together
      character*16 junk
      logical statuspstat
C Common handles
      COMMON/pstathandels/i_mcloop, i_together, junk


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
        ld(1)  = inodesbig  !must be physical dimension of local array
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
     &     mitob(numgroups), 0)
      call ga_pgroup_brdcst(idmedium,msg_dntmc7+MSGINT, groupnums,
     &     mitob(numgroups), 0)

C Destroy GA "proc list"
      if (.not. ga_destroy(g_proc)) call errquit('dntmc_setup_ga:
     &ga_destroy(g_proc) failed', GA_ERR)

C Is done by the new version of util_sgstart
C Create Zero Group
C      idzero = ga_pgroup_create(procnums, util_sgroup_numgroups())

C Create Zeros GAs
      if (ga_pgroup_nodeid(idmedium).eq.0) then
        call ga_pgroup_set_default(idzero)
C allocate timers
      if(util_print('Timers',print_debug)) then
      statuspstat =
     &  pstat_allocate('dntmc:together',Pstat_Qstat,0,junk,i_together)
      endif
C End timers
          ndim = 2
C Changed to dimension (nob,nspc,numgroups)
          dims(1) = NOB*nspc
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
C     &idmedium, idzero, g_pr,  groupnums, numgroups, nob, nspc)
C***********************************************************************
C***********************************************************************
      subroutine dntmc_setup_group
     &(rtdb, idmedium, myid, groupnums, prcfil, prcfil2, fprefix,
     & numgroups)
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
      integer numgroups ! input
      integer groupnums(numgroups) ! input
      character*4 prcfil, prcfil2
      integer i
C Indect output of file units
      character*256 fprefix, fprefixcat ! internal use only !fprefix output
      character*256 fprefix2, dntmc_dir

C Begin Main Program
      call dntmc_write_prcfil(idmedium, groupnums, prcfil, numgroups)
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
      subroutine dntmc_write_prcfil(idmedium,groupnums,prcfil,numgroups)
C***********************************************************************
      Implicit none

C Include Statements
#include "gibbsP.fh"
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
      subroutine dntmc_setup_dntmc(rtdb, iacceptt, iacceptr, iacceptv,
     &  ivtry, ittry, irtry, iacceptpr, iacceptpi, cont1, cont2, nmols,
     &  mprnt, temp, rconf, lmin, lmax, myid, idmedium, groupnums,
     &  prcfil, uu, uavr, usqrd, rsimsum, rsimavr, nspc, dipole,
     &  conv_perc, numgroups, nattotal, natmax, natom,
     &  rinit, uspc, m, spc, tag, mtot, r, rcm, rs, ri, nob, number,
     *  numbersum, numberavr, spctag)
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
      integer mcstps, rsmcstps, mprnt ! mcstps and rsmcstps is internal
      double precision temp
      logical rconf
      integer lmin, lmax
      integer myid, idmedium ! input
      integer node ! internal use only
      integer i, j ! dummy indexes
      double precision rx, ry, rz, dist ! internal use only
      integer molmax, IR, natot ! internal use only
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
      integer numgroups
      integer groupnums(numgroups) ! input

C Variables for Common Blocks
      double precision e, k, esqrd, kb, pi
      double precision disp, alpha
      INTEGER ROTF
      integer idum
      INTEGER AB
      double precision rmin, rmax, rsim
      double precision rclcm
      integer nmin, ncurr, nmax
      double precision dntSEED

C Common Blocks
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/TRPS/DISP,ALPHA,ROTF
      COMMON/dntSEED/IDUM
      COMMON/SEEDS/dntSEED(1)
      COMMON/AORB/AB
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/COORD2/RCLCM(3)
      COMMON/GRAND/NMIN,NCURR,NMAX

C Variable to translate into arrays
      integer nmols, natmax, nattotal, nspc, nob

      double precision rinit(nmols,natmax,3)
      double precision distcalc(nmols) ! internal use only
      character*16 tag(nattotal)
      integer natom(nspc)
      integer uspc(nspc)
      double precision m(nspc,natmax), mtot(nspc)
      integer spc(nmols)
      double precision r(nmols,natmax,3), rcm(nmols,3)
      double precision rs(nmols,natmax,3), ri(nmols,natmax,3)
C Changed dimension to (nob,nspc)
      integer number(nob,nspc)
      double precision numbersum(nob,nspc)
      double precision numberavr(nob,nspc)
      integer spctag(nspc)

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
      dntSEED(1) = 962525439.0d+00

C Aquire Input Information
      if (.not. rtdb_get(rtdb, 'dntmc:norestart', MT_INT, 1, CONT1))
     &   CONT1 = 1

      call util_get_rtdb_state(rtdb, os, oc, or)
      if (or .eqv. .true.) then
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

C Already done in dntmc_allocate_Arrays
C      call dntmc_getspecinfo(rtdb, nmols, natom, rinit, uspc, nspc,
C     &      m, SPC, tag)

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

C Mass of the system
C Initialize MTOT
      Do I = 1, NSPC
        MTOT(USPC(I)) = 0.00D+00
      EnDDo
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

      IF(ga_nodeid() .eq. 0) then
      if(util_print('debug',print_debug)) then
      DO I=1, nspc
        write(6,*)'mtot for species ',i,' is: ',MTOT(I)
        call flush(6)
        DO J = 1, natom(uspc(I))
           write(6,*)'M for species ',i,' and atom ',j,' is: '
     & ,M(USPC(I),J)
           call flush(6)
        EndDo
      Enddo
      endif
      endif
      

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
      Call dntmc_initialstate2(nmols, rinit, radiuscube, natmax,
     & nspc, r, rcm, rs, natom, spc, m, nob, mtot)
C      NODE = idmedium
       NODE = util_sgroup_mygroup()
C Center on COM
      CALL dntmc_CLUSTERCMASS(NMOLS,NMOLS,nspc,natmax,m
     &  ,r,rcm,rs,natom,spc,mtot)
      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO
C find outermost molecule and vector
      DO I = 1, NMOLS
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
C
C Recenter
      CALL dntmc_CLUSTERCMASS(NMOLS,NMOLS,nspc,natmax,m,
     &  r,rcm,rs,natom,spc,mtot)
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
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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

C set spctag
      Do I = 1, nspc
        spctag(I) = 0
C Array of nspc elements.  Labels a configuration
C Each molucule of a particular species adds one to the
C Corresponding element which refers to the species.
      EndDo
      Do I = 1, nmols
        if (DABS(DISTcalc(I) - DIST) .le. 1.0D-10) then
          spctag(spc(I)) = spctag(spc(I)) + 1
        endif
      EndDo

      else
        Do I = 1, nspc
          spctag(I) = 0
        ENDDO
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
        Call dntmc_readrestart(rtdb, uu, uavr, usqrd, dipole,
     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, ivtry, ittry, irtry,
     &  natmax, nob, r, rcm, rs, number, numbersum, numberavr, nspc,
     &  nattotal, natom, spc, tag, spctag)
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
      CALL dntmc_CLUSTERCMASS(NMOLS,NMOLS,nspc,natmax,m,
     &  r,rcm,rs,natom,spc,mtot)
      DO I = 1, NCURR
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO
      DO I = 1, NCURR
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
      call dntmc_setcoords(rtdb, r, nmols, natom, spc, nspc, natmax)
C Get initial energy if first run or new markov chain
      IF((CONT2 .eq. 0).or.(new_markov_chain)) then
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 'Energy did not converge'
          endif
          uu = -0.5D+00                    ! set the energy arbitrarily high - H atom
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
C     &  ivtry, ittry, irtry, iacceptpr, iacceptpi, cont1, cont2, nmols,
C     &  mprnt, temp, rconf, lmin, lmax, myid, idmedium, groupnums,
C     &  prcfil, uu, uavr, usqrd, rsimsum, rsimavr, nspc, dipole,
C     &  conv_perc, numgroups, nattotal, natmax, natom,
C     &  rinit, uspc, m, spc, tag, mtot, r, rcm, rs, ri, spctag)
C***********************************************************************
C***********************************************************************
      Subroutine dntmc_initialstate2(nmols, rinit, radiuscube, natmax,
     & nspc, r, rcm, rs, natom, spc, m, nob, mtot)
C***********************************************************************
      Implicit none
C Include Files
#include "global.fh"
#include "tcgmsg.fh"
#include "gibbsP.fh"
#include "util.fh"

C Variable Declarations
      double precision radiuscube ! input
      integer i, j, k, P
      double precision DELTA
      dimension delta(3)
      double precision dntmc_ran3
      logical overlapflag
      double precision beta, sdbeta, cdbeta
      integer ia1, ia2, ia3
      double precision comp2, comp3, comp2r, comp3r
      integer IR
      double precision dr

      integer nmols ! input
      integer natmax ! input
      integer nspc ! input
      double precision rinit(nmols,natmax,3) !input
      double precision Maxradius(nmols), radius
      double precision r(nmols,natmax,3), rcm(nmols,3)
      double precision rs(nmols,natmax,3)
      integer natom(nspc)
      integer spc(nmols)

      integer nob
      double precision m(nspc,natmax)
      double precision mtot(nspc)

C Common Block Variables
      integer IDUM
      integer ROTF
      double precision RCLCM
      double precision e,ki,esqrd,kb,pi

C Common Blocks
      COMMON/dntSEED/IDUM
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
      CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
      Do J = 1, NATOM(SPC(I))
        Do K = 1,3
          R(I,J,K) = R(I,J,K) - RCM(I,K)
        EndDo
      ENDDO

C Define non-overlapping sphere + 10%
      CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
        CALL dntmc_CMASS(J,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
      CALL dntmc_CLUSTERCMASS(I,NMOLS,nspc,natmax,m,r,rcm,rs,natom,spc
     &  ,mtot)
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
      CALL dntmc_CLUSTERCMASS(I,NMOLS,nspc,natmax,m,r,rcm,rs,natom,spc
     &  ,mtot)
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
        CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
C End of Subroutine dntmc_initialstate2(nmols, rinit, radiuscube, natmax,
C     & nspc, r, rcm, rs, natom, spc, m, nob, mtot)
C***********************************************************************
C***********************************************************************
      SUBROUTINE dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc
     &   ,mtot)
C***********************************************************************
      IMPLICIT none
#include "gibbsP.fh"

      integer i  ! input

      double precision e, k, esqrd, kb, pi

      integer j

      integer nspc, natmax, nmols
      double precision M(nspc,natmax),mtot(nspc)
      double precision R(nmols,natmax,3),rcm(nmols,3)
      double precision rs(nmols,natmax,3)
      integer natom(nspc)
      integer spc(nmols)
 
      COMMON/CONSTANTS/e,k,esqrd,kb,PI

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
C End of SUBROUTINE dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc
C     &     ,mtot)
C***********************************************************************
C***********************************************************************
      SUBROUTINE dntmc_CLUSTERCMASS(NCURR,NMOLS,nspc,natmax,m,r,rcm,rs
     & ,natom,spc,mtot)
C***********************************************************************
      IMPLICIT none
#include "gibbsP.fh"

      double precision rclcm
      double precision MCLTOT
      integer i
      integer ncurr

     
      integer nspc, natmax, nmols
      double precision M(nspc,natmax),mtot(nspc)
      double precision R(nmols,natmax,3),rcm(nmols,3)
      double precision rs(nmols,natmax,3)
      integer natom(nspc)
      integer spc(nmols)

      COMMON/COORD2/RCLCM(3)


      MCLTOT = 0.D+00
      RCLCM(1) = 0.D+00
      RCLCM(2) = 0.D+00
      RCLCM(3) = 0.D+00

      DO I = 1, Ncurr
         MCLTOT = MCLTOT + MTOT(SPC(I))

         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)

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
C End of SUBROUTINE dntmc_CLUSTERCMASS(NCURR,NMOLS,nspc,natmax,m,r,rcm,rs
C      & ,natom ,spc, mtot)
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
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.d-9)
      dimension ma(55)
      data iff /0/
      save iff, mj, ma, mk, inext, inextp
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=iabs(mseed-iabs(idum))
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
      Subroutine dntmc_readrestart(rtdb, uu, uavr, usqrd, dipole,
     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, ivtry, ittry, irtry,
     &  natmax, nob, r, rcm, rs, number, numbersum, numberavr, nspc,
     &  nattotal, natom, spc, tag, spctag)
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
      double precision rmax, rmin
      integer iacceptt, iacceptv, iacceptr, iacceptpr, iacceptpi
      integer ivtry, ittry, irtry
      integer rsmcstps
      double precision rsh
      integer rtdb ! input
      character*4 prcfil !input
      double precision dipole(3) !input

      integer nmols, natmax, nob, nspc, nattotal
      double precision r(nmols,natmax,3), rcm(nmols,3)
      double precision rs(nmols,natmax,3)
C Changed dimensions to (nob,nspc)
      integer number(nob,nspc)
      double precision numbersum(nob,nspc),numberavr(nob,nspc)
C End of Change
      integer natom(nspc), spc(nmols)
      character*16 tag(nattotal)
      integer spctag(nspc)
      double precision dntSEED

C Common Blocks

      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/AVR/USUM,USQRDSUM
      COMMON/SEEDS/dntSEED(1)

C Main Program
      if (.not. rtdb_get(rtdb, 'dntmc:seed:'//prcfil(1:4), MT_DBL, 1,
     &   dntSEED(1)))
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

       call dntmc_read_dntmc_nwchem(rtdb, nmols, prcfil, natmax,
     &  nspc, nattotal, r, rcm, rs, natom, spc, tag)

C      call dntmc_read_dntmc_nwchem(rtdb, nmols, prcfil)

C Changed dimensions to (nob,nspc)
      if (.not. rtdb_get(rtdb, 'dntmc:number:'//prcfil(1:4), MT_INT,
     &    NOB*nspc,number))
     &   call errquit('dntmc_readrestart:  rtdb error for number',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:numbersum:'//prcfil(1:4), MT_DBL,
     &   NOB*nspc,
     &   numbersum))
     &   call errquit('dntmc_readrestart:  rtdb error for numbersum',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:numberavr:'//prcfil(1:4), MT_DBL,
     &   NOB*nspc,
     &   numberavr))
     &   call errquit('dntmc_readrestart:  rtdb error for numberavr',0,
     &   RTDB_ERR)
C End of Change

      if (.not. rtdb_get(rtdb, 'dntmc:rsh:'//prcfil(1:4), MT_DBL, 1, rsh
     &   ))call errquit('dntmc_readrestart:  rtdb error for rsh',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:rsmcstps:'//prcfil(1:4), MT_INT, 1
     &   ,rsmcstps))
     &   call errquit('dntmc_readrestart:  rtdb error for rsmcstps',0,
     &   RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'dntmc:spctag:'//prcfil(1:4), MT_INT,
     &    nspc,spctag))
     &   call errquit('dntmc_readrestart:  rtdb error for spctag',0,
     &   RTDB_ERR)

C End Subroutine
      return
      end
C***********************************************************************
C End of Subroutine dntmc_readrestart(rtdb, uu, uavr, usqrd, dipole,
C     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
C     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, ivtry, ittry, irtry,
C     &  natmax, nob, r, rcm, rs, number, numbersum, numberavr, spctag)
C***********************************************************************
C***********************************************************************
      Subroutine dntmc_read_dntmc_nwchem(rtdb, nmols, prcfil, natmax,
     &  nspc, nattotal, r, rcm, rs, natom, spc, tag)
C***********************************************************************
      Implicit none
C Include Statements
#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "rtdb.fh"

C Variable Declarations
C All Internal
      integer ncent, i, j, k, rdim ! nmols input
      integer rtdb ! input
      character*4 prcfil ! input

      integer nmols, natmax, nspc, nattotal
      double precision  coords(3,nattotal)
      double precision r(nmols,natmax,3), rcm(nmols,3)
      double precision rs(nmols,natmax,3)
      integer natom(nspc)
      integer spc(nmols)
      character*16 tag(nattotal)

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
C End of Subroutine dntmc_read_dntmc_nwchem(rtdb, nmols, prcfil, natmax,
C     &  nspc, nattotal, r, rcm, rs, natom, spc, tag)
C***********************************************************************
C***********************************************************************
      Subroutine dntmc_mcloop(rtdb,temp,UU,nmols,nspc,iacceptt,iacceptr
     &  ,iacceptv,iacceptpr,iacceptpi,uavr,rsimavr,rconf,prcfil,fprefix
     &  ,mprnt,usqrd,rsimsum,idmedium,idzero,cont1,g_pr,rest,flag
     &  ,nda_old,conv_perc,ivtry,ittry,irtry,dipole,idbig,nob,numgroups,
     & natmax,rb,testme,numbersum,numberavr,nda,ndsig,m,mtot,uspc,
     & nattotal,r,rcm,rs,natom,spc,tag,number,spctag)
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
C Timers
#include "pstat.fh"

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
      integer i,j,s,t, z, kp
      double precision conv_perc ! input
      double precision dipole(3) !input/output
      integer movetry
      double precision sumx
      double precision ERRX
      double precision press
      logical inflex
      double precision errrate(nspc), intpress(nspc)
      double precision ratef(nspc), rbf(nspc)
      integer point(nspc)
      double precision acceptratio
      logical dynamic

      integer NOB
      double precision flux(NOB), flux2(NOB,nspc)
      double precision sumerr(NOB), sumerr2(NOB,nspc)
      double precision rflux(NOB,nspc), ratec(NOB,nspc)

      integer numgroups, natmax
C Changed dimension to (nob,nspc)
      double precision nda_old(NOB,nspc) ! input/output
C Changed to dimension (nob,nspc,numgroups)
      double precision rb(nob), testme(nob,nspc,numgroups)
C Changed dimension to (nob,nspc)
      double precision numbersum(nob,nspc), numberavr(nob,nspc)
      integer number(nob,nspc)
C End of change
C Changed dimension to (nob,nspc)
      double precision nda(nob,nspc), ndsig(nob,nspc)
      double precision m(nspc,natmax),mtot(nspc)
      integer uspc(nspc)

      integer nattotal
      double precision r(nmols, natmax, 3), rcm(nmols,3)
      double precision rs(nmols, natmax,3)
      integer natom(nspc)
      integer spc(nmols)
      character*16 tag(nattotal)
      integer spctag(nspc)

      
C Variables for common blocks
      integer L, ncurr,nmin,nmax
      double precision rmin, rmax, rsim
      double precision e,k,esqrd,kb,PI
      double precision disp, alpha
      integer rotf

C  timers common variables
      integer i_mcloop,i_together
      character*16 junk
      logical statuspstat
      integer idbig

C Common Blocks

      COMMON/MSTP/L
      COMMON/GRAND/NMIN,NCURR,NMAX
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/TRPS/DISP,ALPHA,ROTF

C Common handles
      COMMON/pstathandels/i_mcloop, i_together, junk


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
         CALL dntmc_MOVE_RCONF(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,
     +               ACCEPTPR,ACCEPTPI,UU, rtdb, nmols,MOVETRY,dipole,
     & natmax, nspc, nattotal, r, rcm, rs, natom, spc, tag, m, mtot
     &, spctag)
      else
         Do i = 1, nspc
           spctag(I) = 0
         ENDDO
         CALL dntmc_MOVE(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
     +             ,ACCEPTPI,UU, rtdb, nmols, movetry,
     & natmax, nspc, nattotal, r, rcm, rs, natom, spc, tag, m, mtot)
      endif
C         CALL dntmc_ROOTMEANSQUARE(L,NCURR,NSPC)
C Removed ! Became Disconneted

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
         CALL dntmc_NUMBERHIST(L,rconf,nob,number,nspc,spctag)

       RINC = (RMAX - RMIN)/DBLE(NOB)
       RSH = RMIN
       DO I = 1, NOB
          RSH = RSH + RINC
          RB(I) = RSH
       ENDDO

C Calculate Radial Distribution
         CALL dntmc_NUMBERAVERAGE(L,temp,rconf,nob,numgroups,number,
     & numbersum,numberavr,rb,testme,nspc)
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
C Stop Timers
       if(util_print('Timers',print_debug)) then
        Call pstat_off(i_mcloop)
       endif
C End Timers
C Grouping Together Group Information
       if (ga_nodeid().eq.0) then
         call ga_pgroup_set_default(idzero)
C Start Timers
         if(util_print('Timers',print_debug)) then
           call pstat_on(i_together)
         endif
C End Timers

         ilo = 1
C Changed dimension to (nob,nspc)
         ihi = NOB*nspc
         jlo = ga_nodeid() + 1
         jhi = ga_nodeid() + 1

C Changed dimension to (nob,nspc)
         ld(1)  = NOB*nspc   !must be physical dimension of local array
         ld(2)  = 1
         call ga_put(g_pr,ilo,ihi,jlo,jhi,NUMBERAVR,ld)
         call ga_sync
         jlo = 1
         jhi = ga_nnodes()
C Changed to dimension (nob,nspc,numgroups)
         call ga_get(g_pr,ilo,ihi,jlo,jhi,TESTME,ld)
         CALL dntmc_AVERAGEDIST(ga_nnodes(),PERCENT,nda_old,nob,
     &         numgroups,rb, testme, nda, ndsig, nspc)
C Stop Timers
         if(util_print('Timers',print_debug)) then
            call pstat_off(i_together)
         endif
C End Timers
       endif  ! idzero section
       call ga_pgroup_set_default(idmedium)
       call ga_brdcst(msg_dntmc0+MSGDBL, PERCENT, mdtob(1), 0)
       call ga_brdcst(msg_dntmc1+MSGDBL, NDA, mdtob(NOB*nspc), 0)
       call ga_brdcst(msg_dntmc2+MSGDBL, NDSIG, mdtob(NOB*nspc), 0)

       if(ga_pgroup_nodeid(idzero) .eq. 0) then
        if(util_print('fdist',print_low)) then
         WRITE(3,*)(L*util_sgroup_numgroups()), 'Total Configurations'
        endif
        if(util_print('fdist',print_low)) then
         DO Kp = 1, nspc
          WRITE(3,*)'Species number ',Kp
          DO I=1, NOB
           WRITE(3,*)RB(I),NDA(I,Kp),NDSIG(I,Kp)
          ENDDO
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
       Do j = 1, nspc
        Do i=1, NOB
C Changed dimension to (nob,nspc)
         nda_old(i,j) = NDA(i,j)
        enddo
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

C  Rconfig Data Analysis
C  nspcs > 1 not suppported for rsim methodology (results are given but not correct)
       IF (rconf) then
C  Calculate Flux
        flux(1) = 0.00D+00
        Do Kp = 1, nspc
         flux(1) = flux(1)+NDA(1,Kp)
        EndDo
        
        Do i=2, NOB
         flux(i) = flux(i-1)
         Do Kp = 1, nspc
          flux(i) = flux(i) + NDA(i,Kp)
         EndDo
        End Do ! CPF in flux
        Do Kp = 1,nspc
         Do i = 1, NOB
          flux2(i,Kp) = NDA(i,Kp)
         End Do 
        End Do ! Derivative CPF in flux2
C  Errors in Flux
        sumerr(1) = 0.00D+00
        Do Kp = 1, nspc
          sumerr(1) = sumerr(1) + (NDSIG(1,Kp)**2.00D+00)
        EndDo

        Do i = 2, NOB
         sumerr(i) = sumerr(i-1)
         Do Kp = 1, nspc
          sumerr(i) = sumerr(i) + (NDSIG(i,Kp)**2.00D+00)
         EndDo
        END DO
        Do i = 1, NOB
          sumerr(i) = DSQRT(sumerr(i))
        End DO !Error in CPF in sumerr
        Do Kp = 1, nspc
         Do i = 1, NOB
          sumerr2(i,Kp) = NDSIG(i,Kp)
         EndDo
        ENDDO !Error in Derivative CPF in sumerr2
       ELSE
        Do i = 1, nob
          flux(i) = 0.00D+00
          Do Kp = 1, nspc
           flux(i) = flux(i) + NDA(i,Kp)
          EndDo
        EndDo
        Do i = 1, NOB
          flux(i) = flux(i)/(4.D+00*PI*(RB(i)**2.D+00))
          flux(i) = flux(i)*DEXP((4.D+00/3.D+00*PI*(RB(i)**3.D+00))
     &              *press/kb/temp)
        ENDDO ! CPF in flux
        Do Kp = 1, nspc
         Do i = 1, NOB
          if ((i .eq. 1) .or. (i .eq. NOB)) then
            flux2(i,Kp) = 0.00D+00
          else
          flux2(i,Kp) = (flux(i+1) - flux(i-1))/2.00D+00/DBLE(nspc)
          endif
         EndDo
        ENDDO ! Derivative CPF in flux2
C  Errors in Flux
        Do i = 1, NOB
          sumerr(i) = 0.00D+00
          Do Kp = 1, nspc
            sumerr(i) = sumerr(i) + NDSIG(i,Kp)
          ENDDO
        ENDDO
        Do i = 1, NOB
          sumerr(i) = sumerr(i)/(4.D+00*PI*(RB(i)**2.D+00))
          sumerr(i) = sumerr(i)*DEXP((4.D+00/3.D+00*PI*(RB(i)**3.D+00))
     &              *press/kb/temp)
        ENDDO !Error in CPF in sumerr
        Do Kp = 1, nspc
         DO i = 1, NOB
          if ((i .eq. 1) .or. (i .eq. NOB)) then
            sumerr2(i,Kp) = 0.00D+00
          else
           sumerr2(i,Kp)=(sumerr(i-1)**2.00D+00)+(sumerr(i+1)**2.00D+00)
            sumerr2(i,Kp) = DSQRT(sumerr2(i,Kp))/2.00D+00/DBLE(nspc)
          endif
         ENDDO
        ENDDO !Error in Derivative CPF in sumerr2
       ENDIF
C  Calculate ERROR
       ERRX = 0.00D+00
       Do i=1, NOB - 1
        ERRX = ERRX + RB(i+1) - RB(i)
       End Do
       ERRX = ERRX/DBLE(NOB - 1) ! Error related to bin interval
C  Reactive Flux
       Do Kp = 1, nspc
        Do i=1, NOB
         rflux(i,Kp) = flux2(i,Kp)*kb*temp/ERRX
        enddo
       enddo ! Reactive flux in rflux
C Rate Constants
       Do Kp =1, nspc
        Do i=1, NOB
         if(flux(i) .eq. 0.D+00) then
          ratec(i,Kp) = 1.D+99
         else
          ratec(i,Kp) = 1.D+00*rflux(i,Kp)/flux(i)*1.D+13/
     &       DSQRT(2.D+00*PI*kb*temp*mtot(uspc(Kp))/4.184D+00) ! only first species
         endif
        EndDo
       Enddo ! Rate Constants in ratec
C Starting location is Max Flux.  Based on position will go forward or backward
C for minumum search
C Find Minimum Flux
C Minima is point i .lt. i-1 .and. .le. i+1 (i+1 .gt. 0)
C Global minimum is point i .lt. sumx
C Find Maximum Flux
       Do Kp = 1, nspc

        sumx = rflux(1,Kp)
        j = 1
        Do i = 2, NOB
         IF (rflux(i,Kp) .gt. sumx) then
          sumx = rflux(i,Kp)
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
          IF ( (rflux(i,Kp) .lt. rflux(i-1,Kp)) .and.
     &       (rflux(i,Kp) .le. rflux(i+1,Kp)) .and.
     &       (rflux(i,Kp) .lt. sumx)       .and.
     &       (rflux(i+1,Kp) .gt. 0.00D+00) ) then
            sumx = rflux(i,Kp)
            j = i !indux in j
          endif
         End Do

        else

         DO i = z , NOB - 1
          IF ( (rflux(i,Kp) .lt. rflux(i-1,Kp)) .and.
     &       (rflux(i,Kp) .le. rflux(i+1,Kp)) .and.
     &       (rflux(i,Kp) .lt. sumx)       .and.
     &       (rflux(i+1,Kp) .gt. 0.00D+00) ) then
             sumx = rflux(i,Kp)
             j = i ! index in j
          endif
         End Do

        endif

        point(Kp) = j
       EndDo ! Loop over species
C Results
       Do Kp = 1, nspc
        errrate(Kp) = ratec(point(Kp),Kp)*DSQRT((sumerr(point(Kp))
     &          /flux(point(Kp)))**2.D+00
     &          + (flux2(point(Kp),Kp)/sumerr2(point(Kp),Kp))**2.D+00)
        intpress(Kp) = rflux(point(Kp),Kp)/flux(point(Kp))/4/PI
     &          /(RB(point(Kp))**2.00D+00)
        intpress(Kp) = intpress(Kp)*kb*temp
        intpress(Kp)=intpress(Kp)*4.184*1.0D+30/1.01325D+02/6.022142D+23
        ratef(Kp) = ratec(point(Kp),Kp)
        rbf(Kp) = RB(point(Kp))

C Print Results
        if(util_print('Information',print_none)) then
         if(ga_pgroup_nodeid(idzero) .eq. 0) then
          WRITE(6,"(A7,I16, F12.3, ES16.3E3, ES16.3E3, F16.3, F16.3,
     &    F16.3)")
     &    '   @   ',(L*util_sgroup_numgroups()), rbf(Kp), ratef(Kp)
     &    , errrate(Kp)
     &    ,aratiott,aratiotr, percent
          call flush(6)
         endif
        endif
       ENDDO ! End OF nspc loop

       if(util_print('Timers',print_debug)) then
        call ga_pgroup_set_default(idbig)
        call pstat_print_all
        call ga_pgroup_set_default(idmedium)
       endif

C Start Timers
       if (.not. flag) then
        if(util_print('Timers',print_debug)) then
         Call pstat_on(i_mcloop)
        endif
       endif
C End Timers
      ENDIF ! ON MOD(L,mprnt*nob)

C Restart files
      lnow = L + 1
      IF(CONT1.EQ.1)THEN
       if (ga_nodeid() .eq. 0) then
        call dntmc_writefiles(uu, uavr, usqrd, ivtry, ittry, irtry,
     +    rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +    iacceptpi, rsh, prcfil, nmols, lnow, fprefix, dipole, nob,
     &  natmax,nspc, nattotal, r, rcm, rs, number, numbersum, numberavr,
     &  natom, spc, tag, spctag)
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
C     &  ,nda_old,conv_perc,ivtry,ittry,irtry,dipole,idbig,nob,numgroups,
C     & natmax,rb,testme,numbersum,numberavr,nda,ndsig,m,mtot,uspc,
C     & nattotal,r,rcm,rs,natom,spc,tag,number,spctag)
C***********************************************************************
C***********************************************************************
      SUBROUTINE dntmc_MOVE_RCONF(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,
     +               ACCEPTPR,ACCEPTPI,UOLD, rtdb, nmols,MOVETRY,dipold,
     & natmax, nspc, nattotal, r, rcm, rs, natom, spc, tag, m, mtot
     &, spctagold)
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
      integer movetry        ! output
      double precision dipold(3) ! input

      logical task_energy
      external task_energy

      double precision rmin, rsim, rmax
      double precision rclcm
      double precision e, k, esqrd, kb, pi
      double precision disp, alpha
      INTEGER ROTF
      INTEGER NMIN,NCURR,NMAX
      integer l
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

      integer nmols, natmax, nspc, nattotal
      double precision R(nmols,natmax,3),rcm(nmols,3)
      double precision rs(nmols,natmax,3)
      integer natom(nspc)
      integer spc(nmols)
      character*16 tag(nattotal)
      double precision DISTcalc(nmols)

      double precision m(nspc,natmax)
      double precision mtot(nspc)
      integer spctagold(nspc),spctagnew(nspc)

C Added array for configuration files
      integer ti
C End added array
      double precision dntSEED

      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/COORD2/RCLCM(3)
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/TRPS/DISP,ALPHA,ROTF
      COMMON/GRAND/NMIN,NCURR,NMAX
      COMMON/MSTP/L
      COMMON/SEEDS/dntSEED(1)
      COMMON/RNNUMS/RNNUM(1)
      COMMON/dntSEED/IDUM

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

      CALL dntmc_CLUSTERCMASS(NMOLS,NMOLS,nspc,natmax,m
     &  ,r,rcm,rs,natom,spc,mtot)

      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

C  Get new rsim

      DO I = 1, NMOLS
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
C set spctag
      Do I = 1, nspc
        spctagnew(I) = 0
C Array of nspc elements.  Labels a configuration
C Each molucule of a particular species adds one to the
C Corresponding element which refers to the species.
      EndDo
      Do I = 1, nmols
        if (DABS(DISTcalc(I) - DIST) .le. 1.0D-10) then
          spctagnew(spc(I)) = spctagnew(spc(I)) + 1
        endif
      EndDo

C Reject based on volume

C      if ((rnew .gt. rmax).or.(rnew .lt. rmin)) goto 20
       if (rnew .gt. rmax) goto 20

C Do Metropolis on translatonal move

      call dntmc_setcoords(rtdb, r, nmols, natom, spc, nspc, natmax)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 'Energy did not converge'
          endif
          unew = -0.5D+00                    ! set the energy arbitrarily high - H atom
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
      DO I = 1, nspc
        spctagold(I) = spctagnew(I)
      ENDDO

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

      CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)

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

      CALL dntmc_CLUSTERCMASS(NMOLS,NMOLS,nspc,natmax,m
     &  ,r,rcm,rs,natom,spc,mtot)

      DO I = 1, NMOLS
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

C  Get new rsim

      DO I = 1, NMOLS
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
C set spctag
      Do I = 1, nspc
        spctagnew(I) = 0
C Array of nspc elements.  Labels a configuration
C Each molucule of a particular species adds one to the
C Corresponding element which refers to the species.
      EndDo
      Do I = 1, nmols
        if (DABS(DISTcalc(I) - DIST) .le. 1.0D-10) then
          spctagnew(spc(I)) = spctagnew(spc(I)) + 1
        endif
      EndDo

C Reject based on volume

C      if ((rnew .gt. rmax).or.(rnew .lt. rmin)) goto 35
      if (rnew .gt. rmax) goto 35

C Do Metropolis on rotation

      call dntmc_setcoords(rtdb, r, nmols, natom, spc, nspc, natmax)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 'Energy did not converge'
          endif
          unew = -0.5D+00                    ! set the energy arbitrarily high - H atom
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
      DO I = 1, nspc
        spctagold(I) = spctagnew(I)
      ENDDO
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
C End of SUBROUTINE dntmc_MOVE_RCONF(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,
C     +               ACCEPTPR,ACCEPTPI,UOLD, rtdb, nmols,MOVETRY,dipold,
C     & natmax, nspc, nattotal, r, rcm, rs, natom, spc, tag, m, mtot,
C     & spctagold)
C***********************************************************************
C***********************************************************************
      SUBROUTINE dntmc_MOVE(TEMP,ACCEPTT,ACCEPTR,ACCEPTV,ACCEPTPR
     +               ,ACCEPTPI,UOLD, rtdb, nmols, movetry,
     & natmax, nspc, nattotal, r, rcm, rs, natom, spc, tag, m, mtot)
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
      integer movetry        ! output

      logical task_energy
      external task_energy

      double precision rmin, rsim, rmax
      double precision rclcm
      double precision e, k, esqrd, kb, pi
      double precision disp, alpha
      INTEGER ROTF
      INTEGER NMIN,NCURR,NMAX
      integer l
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

      integer nmols, natmax, nspc, nattotal
      double precision R(nmols,natmax,3),rcm(nmols,3)
      double precision rs(nmols,natmax,3)
      integer natom(nspc)
      integer spc(nmols)
      character*16 tag(nattotal)
      double precision DISTcalc(nmols)

      double precision m(nspc,natmax)
      double precision mtot(nspc)

C Added array for configuration files
      integer ti
C End added array
      double precision dntSEED
      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/COORD2/RCLCM(3)
      COMMON/CONSTANTS/e,k,esqrd,kb,PI
      COMMON/TRPS/DISP,ALPHA,ROTF
      COMMON/GRAND/NMIN,NCURR,NMAX
      COMMON/MSTP/L
      COMMON/SEEDS/dntSEED(1)
      COMMON/RNNUMS/RNNUM(1)
      COMMON/dntSEED/IDUM

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
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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

      CALL dntmc_CLUSTERCMASS(NCURR,NMOLS,nspc,natmax,m
     &  ,r,rcm,rs,natom,spc,mtot)

      DO I = 1, NCURR
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

      DO I = 1, NCURR
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
      call dntmc_setcoords(rtdb, r, nmols, natom, spc, nspc, natmax)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 'Energy did not converge'
          endif
          unew = -0.5D+00                    ! set the energy arbitrarily high - H atom
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
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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

      CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)

         DO J = 1, NATOM(SPC(I))
            COMP2 = (R(I,J,IA2) - RCM(I,IA2))
            COMP3 = (R(I,J,IA3) - RCM(I,IA3))

            COMP2R=CDBETA*COMP2+SDBETA*COMP3
            COMP3R=-SDBETA*COMP2+CDBETA*COMP3

            R(I,J,IA2)= RCM(I,IA2) + COMP2R
            R(I,J,IA3)= RCM(I,IA3) + COMP3R
         END DO

      END DO

      CALL dntmc_CLUSTERCMASS(NCURR,NMOLS,nspc,natmax,m
     &  ,r,rcm,rs,natom,spc,mtot)

      DO I = 1, NCURR
         DO J = 1, NATOM(SPC(I))
            R(I,J,1) = R(I,J,1) - RCLCM(1)
            R(I,J,2) = R(I,J,2) - RCLCM(2)
            R(I,J,3) = R(I,J,3) - RCLCM(3)
         END DO
      END DO

      DO I = 1, NCURR
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
      call dntmc_setcoords(rtdb, r, nmols, natom, spc, nspc, natmax)
      if (.not. task_energy(rtdb)) then
        if (.not. task_energy(rtdb)) then   ! try running the energy again
          if(util_print('MCout',print_low)) then
          if(ga_nodeid() .eq. 0) write(12,*) 'Energy did not converge'
          endif
          unew = -0.5D+00                    ! set the energy arbitrarily high - H atom
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
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
         CALL dntmc_CMASS(I,nspc,natmax,nmols,m,r,rcm,rs,natom,spc,mtot)
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
C     +               ,ACCEPTPI,UOLD, rtdb, nmols, movetry,
C     & natmax, nspc, nattotal, r, rcm, rs, natom, spc, tag, m, mtot)
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
***********************************************************************
      SUBROUTINE dntmc_NUMBERHIST(L,rconf,nob,number,nspc,spctag)
C***********************************************************************
      IMPLICIT none
#include "gibbsP.fh"

      integer l  ! input
      logical rconf ! input

      double precision rmin, rsim, rmax
      double precision dntmc_ran3  

      integer  i, j
      double precision INCR

      integer nob, nspc
C Change dimension to (nob, nspc)
      integer number(nob,nspc)
      integer spctag(nspc)
      integer spctagsum, spctagpoint
      integer IDUM

      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/dntSEED/IDUM


      INCR = (RSIM - RMIN)/(RMAX - RMIN)

      IF(L.EQ.1)THEN
         DO I = 1, NOB
           Do J = 1, nspc
C Changed dimension to (nob,nspc)
            NUMBER(I,J) = 0
           EnDDo
         END DO
      endif

      spctagsum = 0
      Do I = 1, nspc
        spctagsum = spctagsum + spctag(I)
        if (spctagsum .eq. 1) spctagpoint = I
      ENDDO

      IF(rconf) then
      I = INT(INCR*DBLE(NOB)) + 1
        if ((nspc .gt. 1).and.(spctagsum .gt. 1)) then
 56       CONTINUE ! redo
          J = INT(DBLE(nspc)*dntmc_ran3(IDUM)) + 1
          if (J .gt. nspc) J = nspc 
          if (spctag(J) .eq. 0) goto 56
        else if (spctagsum .eq. 1) then
          J = spctagpoint
        else
          J = 1
        endif
      else
      I = NINT(INCR*DBLE(NOB))
      J = 1
      endif

      IF (I .le. 0) I = 1
      IF (I .gt. NOB) I = NOB

C Changed dimension to (nob, nspc)
        NUMBER(I,J) = NUMBER(I,J) + 1

      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_NUMBERHIST(L,rconf,nob,number,nspc,spctag)
C***********************************************************************
C***********************************************************************
      SUBROUTINE dntmc_NUMBERAVERAGE(L,temp,rconf,nob,numgroups,number,
     & numbersum,numberavr,rb,testme,nspc)
C***********************************************************************
      IMPLICIT none
#include "gibbsP.fh"
#include "util_sgroup.fh"

      integer l   ! input

      integer  i, j
      double precision e,k,esqrd,kb,PI
      double precision temp, hold
      logical rconf

      integer nob, numgroups, nspc
C Change Dimenstion to (nob,nspc)
      integer number(nob,nspc)
      double precision numbersum(nob,nspc),numberavr(nob,nspc)
      double precision temporary(nob)
C Change Dimension to testme(nob,nspc,numgroups)
      double precision rb(nob), testme(nob,nspc,numgroups)
      
      COMMON/CONSTANTS/e,k,esqrd,kb,PI

      IF(L.EQ.1)THEN
        DO I = 1, NOB
          Do J = 1, nspc
C Change Dimenstion to (nob,nspc)
           NUMBERSUM(I,J) = 0.D+00
           NUMBERAVR(I,J) = 0.D+00
          Enddo
        END DO
      ENDIF

      if (rconf) then
        DO I = 1, NOB
          DO J = 1, nspc
C Change Dimenstion to (nob,nspc)
          NUMBERSUM(I,J) = DBLE(NUMBER(I,J))/DBLE(L)
          NUMBERAVR(I,J) = NUMBERSUM(I,J)
          if (i .gt. 1) then
            NUMBERSUM(I,J) = NUMBERSUM(I,J) + NUMBERSUM(I-1,J)
          endif
          EndDo
        enddo
C Numbersum and Numberavr are not normalized for each species, but for all species
      else
        DO I = 1, NOB
         Do J = 1, nspc
C Change Dimenstion to (nob,nspc)
          NUMBERSUM(I,J) = DBLE(NUMBER(I,J))/DBLE(L)
          NUMBERAVR(I,J) = NUMBERSUM(I,J)
          NUMBERSUM(I,J) = NUMBERSUM(I,J)/(4.D+00*PI*(RB(I)**2.D+00))
         EndDo
        END DO
C Change to find normalization for all species
        Do I = 1, NOB
          temporary(I) = 0.00D+00
        EnDDO
        Do I = 1, NOB
          Do J = 1, nspc
            temporary(I) = temporary(I) + NUMBERSUM(I,J)
          EndDo
        EndDo
        hold = temporary(1)
        hold = NUMBERSUM(1,J)
        DO I = 2, NOB
          if (temporary(I) .ge. hold) then
            hold = temporary(I)
          endif
        EndDo
        DO I = 1, NOB
         Do J = 1, nspc
C Change Dimenstion to (nob,nspc)
          NUMBERSUM(I,J) = NUMBERSUM(I,J)/hold
         End Do
        End Do
      endif

C Numbersum and Numberavr are not normalized for each species, but for all species
      RETURN
      END
C***********************************************************************
C End of SUBROUTINE dntmc_NUMBERAVERAGE(L,temp,rconf,nob,numgroups,number
C     & numbersum,numberavr,rb,testme,nspc)
C***********************************************************************
C***********************************************************************
C  Added NDA_OLD array
      SUBROUTINE dntmc_AVERAGEDIST(NNODES,PERCENT,nda_old,nob,numgroups,
     & rb, testme, nda, ndsig, nspc)
C***********************************************************************
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012

      IMPLICIT none
#include "gibbsP.fh"
#include "util_sgroup.fh"


      integer nnodes ! input parameters
      double precision percent     ! output
      integer i, j, jstar, k

      integer nob, numgroups, nspc
C Changed to dimension (nspc)
      double precision NDSUM(nspc),NDSIGSUM(nspc),NDSGMX(nspc)

C Changed to dimension (nob,nspc)
      double precision nda_old(nob,nspc)
C Changed to dimension (numgroups,nspc)
      double precision meandist(numgroups+1,nspc)
      double precision sumdist(numgroups+1,nspc)

C Changed to dimension (nob,nspc,numgroups)
      double precision rb(nob),testme(nob,nspc,numgroups)
C Changed to dimension (nob,nspc)
      double precision nda(nob,nspc),ndsig(nob,nspc)

C    initialize all values to zero

c  prefer to do NOB
c      DO J = 1, NOB - 1
       DO J = 1, NOB
        Do K = 1, nspc
         NDA(J,K) = 0.D+00
         NDSIG(J,K) = 0.D+00
        EndDo
      ENDDO

       Do K = 1, nspc
        Do J = 1, numgroups+1
         meandist(J,K) = 0.00D+00
         sumdist(J,K) = 0.00D+00
        EndDo
       EndDo

C    find average for each bin across all nodes, for each species

C      DO J = 1, NOB - 1
       Do K = 1, nspc
        DO J = 1, NOB
         NDSUM(K) = 0.D+00
         DO I = 1, NNODES
C Changed to dimension (nob,nspc,numgroups)
            NDSUM(K) = NDSUM(K) + TESTME(J,K,I)
         ENDDO
C Changed to dimension (nob,nspc)
         NDA(J,K) = NDSUM(K)/DBLE(NNODES)
       ENDDO
      ENDDO

C     find sum of distributions
C     find mean of distribution along coordinate
      Do K = 1, nspc
       DO J = 1, NNODES
C Changed to dimension (numgroups,nspc)
        sumdist(J,K) = 0.00D+00
        meandist(J,K) = 0.00D+00
        Do I = 1, NOB
C Changed to dimension (nob,nspc,numgroups)
          sumdist(J,K) = sumdist(J,K) + TESTME(I,K,J)
          meandist(J,K) = meandist(J,K) + (TESTME(I,K,J)*RB(I))
        EndDo
       EndDo
      EndDo
      Do K = 1, nspc
       Do J = 1, NNODES
C Changed to dimension (numgroups,nspc)
        meandist(J,K) = meandist(J,K)/sumdist(J,K)
       EndDo
      EndDo
      Do K = 1, nspc
C Changed to dimension (numgroups,nspc)
       sumdist(numgroups+1,K) = 0.00D+00
       meandist(numgroups+1,K) = 0.00D+00
      Enddo
      Do K = 1, nspc
       if(NNODES .eq. 1) then
        Do I = 1, NOB
C Changed to dimension (numgroups,nspc) or (nob,nspc)
         sumdist(numgroups+1,K) = sumdist(numgroups+1,K) +NDA_OLD(I,K)
         meandist(numgroups+1,K) = meandist(numgroups+1,K)+(NDA_OLD(I,K)
     &     *RB(I))
        EndDo
       else
        Do I = 1, NOB
C Changed to dimension (numgroups,nspc) or (nob,nspc)
         sumdist(numgroups+1,K) = sumdist(numgroups+1,K) + NDA(I,K)
         meandist(numgroups+1,K) = meandist(numgroups+1,K) + (NDA(I,K)
     &     *RB(I))
        EndDo
       endif
      EndDO

      Do K = 1, nspc
C Changed to dimension (numgroups,nspc)
       if(sumdist(numgroups+1,K) .ne. 0.0d+0) then
         meandist(numgroups+1,K) = meandist(numgroups+1,K)/
     &     sumdist(numgroups+1,K)
       endif
      EndDo
C Means saved for each species

C    find variance for each bin across all nodes, and species
C      DO J = 1, NOB - 1
      IF(NNODES .eq. 1) then
       Do K = 1, nspc
        DO J = 1, NOB
C Changed to dimension (nspc) or (nob,nspc)
         NDSIGSUM(K) = 0.D+00 
         NDSIGSUM(K) = NDSIGSUM(K) + (NDA(J,K) - NDA_OLD(J,K))**2.D+00
         if (NDA(J,K) .eq. 0.D+00) then
           NDSIG(J,K) = DSQRT(NDSIGSUM(K)/DBLE(NNODES))
         else
           NDSIG(J,K) = DSQRT(NDSIGSUM(K)/DBLE(NNODES))/NDA(J,K)
C Is now a percent NDSIG
         endif
        ENDDO
       ENDDO
      ELSE
       DO K = 1, nspc
        DO J = 1, NOB
C Changed to dimension (nspc) or (nob,nspc) or (nob,nspc,numgroups)
         NDSIGSUM(K) = 0.D+00
         DO I = 1, NNODES
            NDSIGSUM(K) = NDSIGSUM(K)+(TESTME(J,K,I)-NDA(J,K))**2.D+00
         ENDDO
         if (NDA(J,K) .eq. 0.D+00) then
           NDSIG(J,K) = DSQRT(NDSIGSUM(K)/DBLE(NNODES-1))
         else
           NDSIG(J,K) = DSQRT(NDSIGSUM(K)/DBLE(NNODES-1))/NDA(J,K)
C Is now a percent NDSIG
         endif
        ENDDO
       ENDDO
      ENDIF

C     FInd Variance in the distribution means
      DO K = 1, nspc
C Changed to dimension (nspc) or (nob,nspc) or (numgroups,nspc)
       NDSIGSUM(K) = 0.00D+00
       Do J = 1, NNODES
        NDSIGSUM(K) = NDSIGSUM(K) + (meandist(J,K)
     &     -meandist(numgroups+1,K))**2.D+00
       EndDo
      EndDo
      Do K = 1, nspc
       if (NNODES .gt. 1) then
C Changed to dimension (nspc)
        NDSIGSUM(K) = NDSIGSUM(K)/DBLE(NNODES-1)
       endif
       NDSIGSUM(K) = DSQRT(NDSIGSUM(K))/((RB(2) - RB(1))*DBLE(NOB))
C Is now a percent relative to the width of a bin
      EnDDO

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
C     Is an average over nspc
      PERCENT = 0.00D+00
      Do K = 1, nspc
       PERCENT = PERCENT + NDSIGSUM(K)*100.00D+00
      EndDo
      PERCENT = PERCENT/DBLE(nspc)

      Do K = 1, nspc
       Do J = 1, NOB
C Changed to dimension (nob, nspc)
        NDSIG(J,K) = NDSIG(J,K)*NDA(J,K)
       End Do
      End Do

      RETURN
      END
C***********************************************************************
C  Added NDA_OLD array
C End of SUBROUTINE dntmc_AVERAGEDIST(NNODES,PERCENT,nda_old,nob,numgroups,
C     & rb, testme, nda, ndsig, nspc)
C***********************************************************************
C************************************************************************
      Subroutine dntmc_writefiles(uu, uavr, usqrd, ivtry, ittry, irtry,
     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, fprefix, dipole, nob,
     &  natmax,nspc, nattotal, r, rcm, rs, number, numbersum, numberavr,
     &  natom, spc, tag, spctag)
C***********************************************************************
      Implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "gibbsP.fh"
#include "inp.fh"

      double precision uu, usum, uavr, usqrdsum, usqrd
      double precision rsim, rsimsum, rsimavr
      double precision rmin, rmax, rinc, dipole(3)
      integer iacceptt, iacceptv, iacceptr, iacceptpr, iacceptpi
      integer ivtry, ittry, irtry
      integer rsmcstps, i, j, ir
      double precision rsh
      character*4 prcfil
      character*256 fprefix
      character*256 fprefixcat
      integer ti

      integer nob, nmols, natmax, nspc, nattotal
      double precision r(nmols,natmax,3), rcm(nmols,3)
      double precision rs(nmols,natmax,3)
C Changed to dimension (nob,nspc)
      integer number(nob,nspc)
      double precision numbersum(nob,nspc), numberavr(nob,nspc)
      double precision dntseed

      integer natom(nspc)
      integer spc(nmols)
      character*16 tag(nattotal)
      integer spctag(nspc)

      COMMON/RADIUS/RMIN,RSIM,RMAX
      COMMON/AVR/USUM,USQRDSUM
      COMMON/dntSEEDS/dntSEED(1)


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

         WRITE(20,*)dntSEED(1)
         WRITE(20,*)UU
         WRITE(20,*)USUM
         WRITE(20,*)UAVR
         WRITE(20,*)USQRDSUM
         WRITE(20,*)USQRD
         WRITE(20,*)(DIPOLE(IR), IR=1,3)
         WRITE(20,*)RMIN, RMAX
         WRITE(20,*)RSIM
         WRITE(20,*)(spctag(I),I=1,nspc)
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
C Changed to dimension (nob,nspc)
            WRITE(20,*)(NUMBER(I,J),J=1,nspc)
            WRITE(20,*)(NUMBERSUM(I,J),J=1,nspc)
            WRITE(20,*)RSH,(NUMBERAVR(I,J),J=1,nspc)
         ENDDO

         CLOSE(20)

      return
      end

C***********************************************************************
C  End of Subroutine dntmc_writefiles(uu, uavr, usqrd, ivtry, ittry, irtry,
C     +  rsimsum, rsimavr, iacceptt, iacceptv, iacceptr, iacceptpr,
C     +  iacceptpi, rsh, prcfil, nmols, rsmcstps, fprefix, dipole, nob,
C     &  natmax, nspc, nattotal, r, rcm, rs, number, numbersum, numberavr,
C     &  natom, spc, tag, spctag)
C***********************************************************************

