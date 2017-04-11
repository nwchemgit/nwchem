      logical function hnd_property(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "cfock.fh"
#include "prop.fh"
#include "inp.fh"
c
c     Property module.
c
c     Does energy run to check if energy was converged far enough,
c     then will do the properties. Some checking. If we want to do 
c     spin-spin or shielding we need to make the geometry C1 during
c     the property calculation. After that we have to reset it back
c     to the proper symmetry.
c
c     1. Initialize all necessary SCF/DFT information
c     2. Do energy
c     3. Get geometry and basis
c     4. Run properties
c
c ... jochen: got terrible problems with SCF convergence, therefore
c             I relaxed the tightened criteria here.
c
      logical  int_normalize
      external int_normalize

      integer rtdb              ! [input] database handle
      integer geom, basis       ! handles
      logical status, usec1
      character*255 title, vectors_in
      character*32 theory
      character*32 operation
      character*8 gridtype
      integer giao, spinspin, geomold, noaoresponse
      integer gshift
      double precision ithresh, tthresh, e_conv, tol_rho, tol_rho_rd
      double precision job_grid_acc
c
      character*32 key
      logical flag
c 
      logical grid_clinit, task_energy
      external grid_clinit, task_energy
c
      logical nocheck
      logical lcgmin

c     Initialization
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
c
c     Get theory type
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('hnd_property: theory not specified',555,
     &       INPUT_ERR)
      if (theory.eq.'hyb') theory = 'dft'
      use_theory = theory
c
c     get external movecs if any (MV)
c     ------------------------------
      if (.not. rtdb_cget(rtdb, 'prop:vectors', 1, vectors_in)) then
        vectors_in = ' '
      end if
c
c     skip checking if requested(MV)
c     ------------------------------
      if (.not. rtdb_get(rtdb, 'prop:nocheck',mt_log, 1, nocheck)) then
        nocheck = .false.
      end if

c ... jochen: try nocheck = .true. always (debugging)
c      nocheck = .true.
c ... jochen: we now set prop:nocheck in prop_input.F for response calc.
c      
c
c     Check if we are running spinspin or shielding. In that case we
c     have to remove symmetry from the geometry.
c
      gshift=0 ! FA-09-16-10
      status=rtdb_get(rtdb,'prop:gshift',MT_INT,1,gshift) ! FA-09-16-10     
      giao = 0
      spinspin = 0
      noaoresponse = 0 
      status=rtdb_get(rtdb,'prop:giao',MT_INT,1,giao)             
      status=rtdb_get(rtdb,'prop:spinspin',MT_INT,1,spinspin)     
      status=rtdb_get(rtdb,'prop:aoresponse',MT_INT,1,noaoresponse)
c     usec1 = giao.eq.0 .or. spinspin.eq.0 .or. noaoresponse.eq.0
c     FA-08-18-10
      usec1 = gshift.eq.0 .or. giao.eq.0 .or. 
     &        spinspin.eq.0 .or. noaoresponse.eq.0

      if (.not.rtdb_get(rtdb, 'dft:cgmin', mt_log, 1, lcgmin))
     &    lcgmin = .false.
c
c     Get geometry and make C1 if necessary
c
      if (usec1) then
         if (.not. geom_create(geomold, 'geometry'))
     $      call errquit('hnd_property: geom_create?', 0, GEOM_ERR)
         if (.not. geom_rtdb_load(rtdb, geomold, 'geometry'))
     $      call errquit('hnd_property: no geometry ', 0, RTDB_ERR)
         if (.not. geom_makec1(geomold, geom)) call
     $    errquit('hnd_property: failed decrease symmetry',555,GEOM_ERR)
         if (.not.geom_rtdb_store(rtdb,geom,'geometry')) call
     &    errquit('hnd_property: geom_rtdb_store failed',555,RTDB_ERR)
         if (.not.geom_destroy(geom)) call
     &    errquit('hnd_property: geom_destroy failed',555,RTDB_ERR)
      endif
c
c     Some printout
c
      if (ga_nodeid().eq.0) then
         call util_flush(luout)
         if (util_print('information',print_low)) then
            call util_print_centered(luout,'NWChem Property Module',
     $           40,.true.)
            write(luout,*)
            call util_flush(luout)
          endif
          if (util_print('information',print_medium)) then
            write(luout,*)
            if (title .ne. ' ') then
               call util_print_centered(luout, title, 40, .false.)
               write(luout,*)
            endif
            call util_flush(luout)
         endif
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom)) 
     $           call errquit('property: geom_print ?',0, GEOM_ERR)
            call util_flush(luout)
         endif
         if (util_print('basis',print_high)) then
            if (.not. bas_print(basis))
     $           call errquit('property: bas_print ?',0, BASIS_ERR)
            call util_flush(luout)
         endif
      endif
c
c     Set up SCF and DFT stuff for shielding and spinspin CPHF
c     Also for general more improved accuracy of properties
c
      if (theory.eq.'dft' .and. .not. nocheck) then
c
c     set grid to fine or leave the user to set it to xfine
c
         if (rtdb_cget(rtdb, 'dft:gridtype', 1, gridtype)) then
            if(gridtype.eq.'coarse'.or.gridtype.eq.'xcoarse'.
     $         or.gridtype.eq.'medium') gridtype='fine'
         else
            gridtype='fine'
         endif
c
         if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     $       call errquit('hnd_property:rtdbcput failed',0, RTDB_ERR)
c
c     set grid_acc consistent with fine (at least), ie 1d-7
c
         if (.not. rtdb_get(rtdb,'dft:job_grid_acc', mt_dbl, 1,
     &      job_grid_acc)) then
            job_grid_acc = 1d-7
         else
            job_grid_acc = min(job_grid_acc,1d-7)
         endif
         if (.not. rtdb_put(rtdb,'dft:job_grid_acc', mt_dbl, 1,
     &       job_grid_acc))call errquit('hnd_property:rtdbput failed',
     &       1,RTDB_ERR)
c
c     weight derivatives bit turned off
c
         if (.not. rtdb_put(rtdb,'dft:weight derivatives', mt_log,1,
     &       .false.)) call errquit('hnd_property:rtdbput failed',2,
     &                               RTDB_ERR)
c
c     set tol_rho to at least 1d-11 (it is 1d-10 otherwise)
c
         tol_rho=1d-11
         if (rtdb_get(rtdb, 'dft:tol_rho', mt_dbl, 1, tol_rho_rd))
     .       tol_rho=min(tol_rho,tol_rho_rd)
         if (.not.rtdb_put(rtdb, 'dft:tol_rho', mt_dbl, 1, tol_rho))
     .        call errquit(' hnd_property: rtdbput failed ',11,0)
c
c     kill symmetry for XC grid build when doing C1 symmetry CPHF
c
         if (usec1) then
           if (.not. rtdb_put(rtdb, 'dft:gridsym',
     $         MT_LOG, 1, .false.)) call errquit(
     $         ' hnd_property: rtdb_put failed',0, RTDB_ERR)
         endif
c
c     get xfac for doing xc_gotxc right on restart
c     but ignore the status since there are not xfac
c     on rtdb when we are in start mode
c
         status=grid_clinit(rtdb)
         call grid_cleanup(.false.)
c
c     set deltaE at least to 1d-7
c
         if (.not.rtdb_get(rtdb,'dft:e_conv',mt_dbl,1,e_conv))then
             e_conv=1d-7
         else
             e_conv=min(e_conv,1d-7)
         endif
         if (.not. rtdb_put(rtdb, 'dft:e_conv', mt_dbl, 1, e_conv))
     &       call errquit('hnd_property: rtdb_put failed', 100,
     &                    RTDB_ERR)
c
         call dft_inpanae(rtdb)
      endif  ! nocheck
c
c     Make sure the energy and vectors have been calculated
c
      if(vectors_in.eq.' ') then
c
c       Remove the converged flag to force recalculation
c
        key = theory(1:inp_strlen(theory))//":converged"
        if (rtdb_get(rtdb, key, mt_log, 1, flag)) then
          if (.not. rtdb_put(rtdb, key, mt_log, 1, .false.))
     $      call errquit("task: could not write to rtdb",0,RTDB_ERR)
        endif
c
c       Calculate the energy
c
        if(.not.task_energy(rtdb)) call errquit
     *    ('hnd_property: energy failure', 555, CALC_ERR)
c
c       Depending on calculation write prop:input vectors to rtdb
c
        if (theory.eq.'dft') then
           if (.not. rtdb_cget(rtdb,'dft:output vectors',1,vectors_in))
     &         call errquit('hnd_property: rtdb_cget failed', 100,
     &                      RTDB_ERR)
        elseif (theory.eq.'scf') then
           if (.not. rtdb_cget(rtdb,'scf:output vectors',1,vectors_in))
     &         call errquit('hnd_property: rtdb_cget failed', 100,
     &                      RTDB_ERR)
        endif
      end if

      if (.not. rtdb_cput(rtdb,'prop:input vectors',1,vectors_in))
     &    call errquit('hnd_property: rtdb_cput failed', 100,
     &                 RTDB_ERR)
c
c    Get geometry and basis set
c
      if (.not. geom_create(geom, 'geometry'))
     $    call errquit('hnd_property: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $    call errquit('hnd_property: no geometry ', 0, RTDB_ERR)
c
c     Get the basis set
c
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('property: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('property: no ao basis ', 0, RTDB_ERR)
      if (.not. int_normalize(rtdb,basis))
     $     call errquit('property: int_normalize', 0, INT_ERR)
c
c     call hondo property calculation routines
c
      if (.not.bas_any_gcorsp(basis)) then
          call hnd_propty ( rtdb, basis, geom )
      else
          write(luout,99991)
      endif
c
      if (.not.(
     &    (bas_destroy(basis))
     &    .and.
     &    (geom_destroy(geom))
     &    ))
     &    call errquit
     &    ('hnd_property:error destroying geom and basis handles',911,
     &       GEOM_ERR)
c
c     Reset geometry if necessary (i.e. after C1 run)
c
      if (usec1) then
         if (.not.geom_rtdb_store(rtdb,geomold,'geometry')) call
     &    errquit('hnd_property: geom_rtdb_store failed',555,RTDB_ERR)
         if (.not.geom_destroy(geomold)) call
     &    errquit('hnd_property: geom_destroy failed',555,RTDB_ERR)
      endif
c
      hnd_property = .true.
c
      return
99991 format(' property routines lifted from -hondo- do not handle',
     &       ' sp functions or general contractions ... .',/,
     &       ' please use the "segment" keyword on the basis set ',/,
     &       ' skip and continue ')
      end
