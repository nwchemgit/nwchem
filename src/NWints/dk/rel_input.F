      subroutine rel_input(rtdb)
*
* $Id: rel_input.F,v 1.16 2008-05-13 20:38:04 niri Exp $
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
c     
c     Read input for the relativistic block
c     
      integer rtdb ! [input] handle to database
c     
      character*255 test
      double precision clight
      logical dk_run             ! Douglas-Kroll
      logical dmd_run            ! Dyall modified Dirac
      logical zora_run           ! ZORA

      logical status
      integer irelKey            ! Relativistic method keys

c     Initialise relativistic approach key
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
     $               'calculations is ',F16.10, 'au.'/)
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
      elseif (inp_compare(.false.,'zora', test)) then
        zora_run = .true.
        if (.not. rtdb_put(rtdb, 'zora', mt_log, 1, zora_run))
     &     call errquit('rel_input: rtdb_put failed', 0, RTDB_ERR)
         call zora_input(rtdb)
c
      elseif(inp_compare(.false.,'zora:cutoff',test)) then
         if(.not. inp_f(zoracutoff)) zoracutoff = 1.0d-08  ! set a default cutoff
         if (.not. rtdb_put(rtdb,'zora:cutoff',mt_dbl,1,zoracutoff))
     &    call errquit('rel_input: rtdb put failed',0, RTDB_ERR)
c
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
