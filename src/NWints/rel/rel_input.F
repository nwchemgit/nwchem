      subroutine rel_input(rtdb)
*
* $Id: rel_input.F 21948 2012-02-03 19:27:57Z niri $
*
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "rel_consts.fh"
#include "errquit.fh"
#include "zora.fh"
#include "modelpotential.fh"
#include "x2c.fh"
c     
c     Read input for the relativistic block
c     
      integer rtdb ! [input] handle to database
      integer slc_spinpolAO ! for selecting spin-density matrix contrib
c     
      character*255 test
      double precision clight
      logical dk_run              ! Douglas-Kroll
      logical dmd_run             ! Dyall modified Dirac
      logical zora_run            ! ZORA
      logical x2c_run             ! exact decoupling (X2C)
      integer mpoption
c
      logical status
      integer irelKey            ! Relativistic method keys
      logical do_prntNMRCS
      logical skip_csAOev,skip_gshiftAOev,
     &        skip_hypAOev,skip_efgz4AOev,
     &        Knucl
      logical skip_cphf_ev_shield,
     &        skip_cphf_ev_gshift,
     &        skip_cphf_ev_hyp

      irelKey = 0
c
c     Set the available relativistic approaches to .false.
c
c     Douglas-Kroll 
      dk_run = .false.
      if (.not. rtdb_put(rtdb, 'doug_kroll', mt_log, 1, dk_run))
     $   call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
c
c     Dyall Modified Dirac
      dmd_run = .false.
      if (.not. rtdb_put(rtdb, 'dyall_mod_dir', mt_log, 1, dmd_run))
     $   call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
c
c     ZORA
      zora_run = .false.
      if (.not. rtdb_put(rtdb, 'zora', mt_log, 1, zora_run))
     &  call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
c
c     X2C
      x2c_run = .false.
      if (.not. rtdb_put(rtdb, 'x2c', mt_log, 1, x2c_run))
     &  call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
c
c     Set some relativistic variables we could need later on
c
   1  if (.not. inp_read()) call errquit('rel_input: inp_read failed',0,
     &       INPUT_ERR)
      if (.not. inp_a(test))
     $     call errquit('rel_input: failed to read keyword', 0,
     &       INPUT_ERR)
c
      if (inp_compare(.false.,'clight', test)) then
c
         if (.not. inp_f(clight)) then
            call errquit('rel_input: failed reading clight', 0,
     &       UNKNOWN_ERR)
         else 
           if(.not.rtdb_put(rtdb,'relativistic:clight',mt_dbl,1,clight))
     $        call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
         endif
         write(luout,10) clight
  10        format(/,'The speed of light used in the relativistic ',
     $               'calculations is ',E15.7, 'au.'/)
c
      elseif (inp_compare(.false.,'douglas-kroll', test)) then
         dk_run = .true.
         if (.not. rtdb_put(rtdb, 'doug_kroll', mt_log, 1, dk_run))
     $     call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
         call dk_input(rtdb)
c
      elseif (inp_compare(.false.,'dyall-mod-dirac', test)) then
         dmd_run = .true.
         if (.not. rtdb_put(rtdb, 'dyall_mod_dir', mt_log, 1, dmd_run))
     $     call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
         call dmd_input(rtdb)
c
      elseif (inp_compare(.false.,'x2c', test)) then
        x2c_run = .true.
        if (.not. rtdb_put(rtdb, 'x2c', mt_log, 1, x2c_run))
     &     call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
         call x2c_input(rtdb)
c
      elseif (inp_compare(.false.,'zora', test)) then
         call zora_input(rtdb)
c
      elseif (inp_compare(.false.,'modelpotential', test)) then
        use_modelpotential = .true.
        if (.not. rtdb_put(rtdb, 'modelpotential', mt_log, 1, 
     &  use_modelpotential))
     &     call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
        if (.not. inp_i(mpoption)) mpoption = 1 ! default is modbas.4c
        if (.not. rtdb_put(rtdb,'modelpotential:option', 
     &                      mt_int,1,mpoption))
     &     call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
        call modelpotential_input(rtdb)
c
      elseif(inp_compare(.false.,'zora:cutoff',test)) then
         if(.not. inp_f(zoracutoff)) zoracutoff = 1.0d-15  ! set a default cutoff
         if (.not. rtdb_put(rtdb,'zora:cutoff',mt_dbl,1,zoracutoff))
     &    call errquit('rel_input: rtdb put failed',0, RTDB_ERR)
c
      elseif(inp_compare(.false.,'x2c:cutoff',test)) then
         if(.not. inp_f(x2ccutoff)) x2ccutoff = 1.0d-15  ! set a default cutoff
         if (.not. rtdb_put(rtdb,'x2c:cutoff',mt_dbl,1,x2ccutoff))
     &    call errquit('rel_input: rtdb put failed',0, RTDB_ERR)
c
c     zoracutoff_EFG, used in get_rhoS() from dft_zora_utils.F
      elseif(inp_compare(.false.,'zora:cutoff_EFG',test)) then
         if(.not. inp_f(zoracutoff_EFG)) zoracutoff_EFG = 1.0d-04 ! set dflt cutoff
         if (.not. rtdb_put(rtdb,'zora:cutoff_EFG', 
     &                      mt_dbl,1,zoracutoff_EFG))
     &    call errquit('rel_input: rtdb put failed',0, RTDB_ERR)
c
c     zoracutoff_NMR, used in get_NMR() from dft_zora_utils.F
      elseif(inp_compare(.false.,'zora:cutoff_NMR',test)) then
         if(.not. inp_f(zoracutoff_NMR)) zoracutoff_NMR = 1.0d-08 ! set dflt cutoff
         if (.not. rtdb_put(rtdb,'zora:cutoff_NMR',
     &                      mt_dbl,1,zoracutoff_NMR))
     &    call errquit('rel_input: rtdb put failed',0, RTDB_ERR)
c
      elseif(inp_compare(.false.,'zora:do_NonRel',test)) then
         do_NonRel = .true.                                  ! set do_NonRel
         if (.not. rtdb_put(rtdb,'zora:do_NonRel', mt_log,1,do_NonRel))
     &    call errquit('rel_input: rtdb put failed',0, RTDB_ERR)
c
      elseif(inp_compare(.false.,'zora:do_prntNMRCS',test)) then
         do_prntNMRCS = .true.                            ! set do_prntNMRCS
         if (.not. rtdb_put(rtdb,'zora:do_prntNMRCS',
     &                      mt_log,1,do_prntNMRCS))
     &    call errquit('rel_input: rtdb put failed',0, RTDB_ERR)
c
      elseif(inp_compare(.false.,'zora:not_zora_scale',test)) then 
         not_zora_scale = .true.                               ! set not_zora_scale
         if (.not. rtdb_put(rtdb,'zora:not_zora_scale',
     &                     mt_log,1,not_zora_scale))
     &    call errquit('rel_input: rtdb put failed',0, RTDB_ERR)
c
c ------- SKIP ev. AO matrices for NMRCS or gshifts or hyperfine
c Note.- It will read instead data from a file, if the file does not
c        exist it will stop calc. with a warning message.
      elseif(inp_compare(.false.,'zora:skip_csAOev',test)) then 
         skip_csAOev = .true.      
         if (.not. rtdb_put(rtdb,'zora:skip_csAOev',           
     &                     mt_log,1,skip_csAOev))             
     &    call errquit('rel_input: rtdb put failed skip csAOev',
     &                  0, RTDB_ERR) 
      elseif(inp_compare(.false.,'zora:skip_gshiftAOev',test)) then
         skip_gshiftAOev = .true.                               
         if (.not. rtdb_put(rtdb,'zora:skip_gshiftAOev',
     &                     mt_log,1,skip_gshiftAOev))
     &    call errquit('rel_input: rtdb put failed skip gshiftAOev',
     &                  0, RTDB_ERR)
      elseif(inp_compare(.false.,'zora:skip_hypAOev',test)) then  
         skip_hypAOev = .true.
         if (.not. rtdb_put(rtdb,'zora:skip_hypAOev',
     &                     mt_log,1,skip_hypAOev))
     &    call errquit('rel_input: rtdb put failed skip hypAOev',
     &                  0, RTDB_ERR)
      elseif(inp_compare(.false.,'zora:skip_efgz4AOev',test)) then   
         skip_efgz4AOev = .true.
         if (.not. rtdb_put(rtdb,'zora:skip_efgz4AOev',
     &                     mt_log,1,skip_efgz4AOev))
     &    call errquit('rel_input: rtdb put failed skip efgz4AOev',
     &                  0, RTDB_ERR)
c
c     definitions for skip-cphf flags
      elseif(inp_compare(.false.,'zora:skip_cphf_ev_shield',test)) then
         skip_cphf_ev_shield = .true.
         if (.not. rtdb_put(rtdb,'zora:skip_cphf_ev_shield',
     &                     mt_log,1,skip_cphf_ev_shield))
     &    call errquit('rel_input-1: rtdb put failed skip cphf shield ',
     &                  0, RTDB_ERR)
      elseif(inp_compare(.false.,'zora:skip_cphf_ev_gshift',test)) then
         skip_cphf_ev_gshift = .true.
         if (.not. rtdb_put(rtdb,'zora:skip_cphf_ev_gshift',
     &                     mt_log,1,skip_cphf_ev_gshift))
     &    call errquit('rel_input-2: rtdb put failed skip cphf gshift',
     &                  0, RTDB_ERR)
      elseif(inp_compare(.false.,'zora:skip_cphf_ev_hyp',test)) then
         skip_cphf_ev_hyp = .true.
         if (.not. rtdb_put(rtdb,'zora:skip_cphf_ev_hyp',
     &                     mt_log,1,skip_cphf_ev_hyp))
     &    call errquit('rel_input-3: rtdb put failed skip cphf hyp',
     &                  0, RTDB_ERR)
c
      elseif(inp_compare(.false.,'zora:Knucl',test)) then
c      Knucle purpose: For including ONLY nuclear part in evaluation of K ZORA
         Knucl = .true. 
         if (.not. rtdb_put(rtdb,'zora:Knucl',
     &                     mt_log,1,Knucl))
     &    call errquit('rel_input: rtdb put failed Knucl',
     &                  0, RTDB_ERR)
      elseif(inp_compare(.false.,'zora:slc_spinpolAO',test)) then
c   For selecting spin-density matrix (A-B means total spin density matrix)
c   0=A-B,1=A,2=-B selecting spin-density matrix  
         if (.not. inp_i(slc_spinpolAO)) slc_spinpolAO= 0  ! set a default 
         if (.not. rtdb_put(rtdb,'zora:slc_spinpolAO',
     &                      mt_int,1,slc_spinpolAO))
     &    call errquit('rel_input: rtdb put failed',0, RTDB_ERR)
      else if (inp_compare(.false.,'end', test)) then
c
c        end will be catched in while loop structure
c
      else
         call errquit('rel_input: unrecognized directive', 0, INPUT_ERR)
      endif
c
c     End of while loop over the input
c
      if (.not. inp_compare(.false.,'end', test)) goto 1
c
c     User could have set multiple relativistic approaches in the input file by accident. 
c     Check if this is the case, stop if multiple methods are set to true.
c
      irelKey = 0
      status = rtdb_get(rtdb, 'doug_kroll', mt_log, 1, dk_run)
      if (dk_run) irelKey = irelKey + 1
c
      status = rtdb_get(rtdb, 'dyall_mod_dir', mt_log, 1, dmd_run)
      if (dmd_run) irelKey = irelKey + 1
c
      status = rtdb_get(rtdb, 'zora', mt_log, 1, zora_run)
      if (zora_run) irelKey = irelKey + 1
c
c     Check if multiple approaches have been set 
      if (irelKey .gt. 1) then
         if (ga_nodeid() .eq. 0) then
           write(luout,20)
  20       format(/'Input error relativistic:'/,
     & 'Multiple relativistic methods have been specified 
     & in the input!')
         endif
         call errquit('rel_input: input error',0, INPUT_ERR)
      endif

      return
      end
