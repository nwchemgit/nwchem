      subroutine gw_input(rtdb)
c
c     Input parser for GW@DFT
c
c     Daniel Mejia-Rodriguez
c     November 2020
c
      implicit none
c
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      integer rtdb
      character*255 test
c
      integer, parameter :: num_dirs = 13
      character*11 dirs(num_dirs)
      data dirs /'evgw0', 'evgw' , 'method', 'solver',
     &           'convergence', 'states', 'first', 'ncap', 
     &           'eta', 'rpa', 'diag', 'core', 'end'/
c
      integer ind, qpiter,eviter,noqp(2),nvqp(2), ngl, first
      double precision threshold, eta
      double precision, parameter :: ha2ev = 27.211386245988d0
c
c     Set default values and put them in RTDB
c
      call gw_defaults(rtdb)
c
c     Parse the input file
c
 10   if (.not. inp_read()) then
        call errquit('gw_input: failed to read input',0,INPUT_ERR)
      endif
c
      if (.not. inp_a(test)) then
        call errquit('gw_input: failed to read keyword',0,INPUT_ERR)
      endif
c
      if (.not.inp_match(num_dirs, .false., test, dirs, ind)) then
        call errquit('gw_input: unknown directive', 0, INPUT_ERR)
      end if
c
      goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 
     $      1100, 1200, 9999) ind

      call errquit('gw_input: unknown directive',ind, INPUT_ERR)
c
c     -----------------------------------------------------------------
c
c     evGW0
c
 100  if (.not. inp_i(eviter)) then
        write(LuOut,*) ' gw_input: evGW0 iterations value not found;',
     &                 ' defaulting to 4'
        eviter = 4
      endif
      if (.not.rtdb_put(rtdb, 'gw:evgw0', mt_log, 1, .true.))
     $  call errquit('gw_input: rtdb_put failed', 100, RTDB_ERR)
      if (.not.rtdb_put(rtdb, 'gw:eviter', mt_int, 1, eviter))
     &  call errquit('gw_input: rtdb_put failed', 100, RTDB_ERR)
      goto 10
c
c     evGW
c
 200  if (.not. inp_i(eviter)) then
        write(LuOut,*) ' gw_input: evGW iterations value not found;',
     &                 ' defaulting to 4'
        eviter = 4
      endif
      if (.not.rtdb_put(rtdb, 'gw:evgw', mt_log, 1, .true.))
     $  call errquit('gw_input: rtdb_put failed', 100, RTDB_ERR)
      if (.not.rtdb_put(rtdb, 'gw:eviter', mt_int, 1, eviter))
     &  call errquit('gw_input: rtdb_put failed', 100, RTDB_ERR)
      goto 10
c
c     method
c
 300  if (.not. inp_a(test))
     &  call errquit('gw_input: method not found', 300, INPUT_ERR)
      if (inp_compare(.false.,'cdgw',test)) then
        if (.not.rtdb_put(rtdb, 'gw:cdgw', mt_log, 1, .true.))
     &    call errquit('gw_input: rtdb_put failed', 400, RTDB_ERR)
        if (.not.inp_i(ngl)) goto 10
        if (.not.rtdb_put(rtdb, 'gw:ngl', mt_int, 1, ngl))
     &    call errquit('gw_input: rtdb_put failed', 300, RTDB_ERR)
!
!      Analytic continuation method not yet implemented
!
!      elseif (inp_compare(.false.,'acgw',test)) then
!        if (.not.rtdb_put(rtdb, 'gw:acgw', mt_log, 1, .true.))
!     &    call errquit('gw_input: rtdb_put failed', 400, RTDB_ERR)
!        call errquit('gw_input: AC-GW not ready yet',0,INPUT_ERR)
!
      elseif (.not. inp_compare(.false.,'analytic',test)) then
        call errquit('gw_input: unknown GW method',300,INPUT_ERR)
      endif
      goto 10
c
c     solver
c
 400  if (.not. inp_a(test)) goto 8000
      if (inp_compare(.false.,'graph',test)) then
        if (.not.rtdb_put(rtdb, 'gw:graph', mt_log, 1, .true.))
     &    call errquit('gw_input: rtdb_put failed', 400, RTDB_ERR)
      elseif (inp_compare(.false.,'newton',test)) then
        if (.not. inp_i(qpiter)) goto 8000
        if (.not.rtdb_put(rtdb, 'gw:newton', mt_log, 1, .true.))
     &    call errquit('gw_input: rtdb_put failed', 400, RTDB_ERR)
        if (.not.rtdb_put(rtdb, 'gw:qpiter', mt_int, 1, qpiter))
     &    call errquit('gw_input: rtdb_put failed', 400, RTDB_ERR)
      else
        goto 8000
      endif
      goto 10
c
c     convergence
c
 500  if (.not. inp_f(threshold))
     &  call errquit('convergence <float> [ev || hartree]',0,INPUT_ERR)
      if (inp_a(test)) then
        if (inp_compare(.false.,'ev',test)) then
          threshold = threshold/ha2ev
        elseif (.not. inp_compare(.false.,'hartree',test)) then
          call errquit('gw_input: invalid energy unit',0,INPUT_ERR)
        endif
      else
        write(LuOut,*) ' gw_input: No convergence unit was specified;',
     &                 ' assuming eV'
      endif
      if (.not. rtdb_put(rtdb, 'gw:threshold', mt_dbl, 1, threshold))
     &  call errquit('gw_input: rtdb_put failed', 500, RTDB_ERR)
      goto 10
c
c     states
c
 600  if (.not. inp_a(test)) goto 8200

      if (inp_compare(.false.,'alpha',test)) then
        if (.not.inp_a(test)) goto 8200
 601    if (inp_compare(.false.,'occ',test)) then
          if (.not. inp_i(noqp)) goto 8200
          if (.not.rtdb_put(rtdb, 'gw:noqp_alpha',mt_int,1,noqp))
     $      call errquit('gw_input: rtdb_put failed', 600, RTDB_ERR)
        elseif (inp_compare(.false.,'vir',test)) then
          if (.not. inp_i(nvqp)) goto 8200
          if (.not.rtdb_put(rtdb, 'gw:nvqp_alpha',mt_int,1,nvqp))
     $      call errquit('gw_input: rtdb_put failed', 600, RTDB_ERR)
        else
          goto 8200
        endif
        if (.not. inp_a(test)) goto 10
        goto 601

      elseif (inp_compare(.false.,'beta',test)) then
        if (.not.inp_a(test)) goto 8200
 602    if (inp_compare(.false.,'occ',test)) then
          if (.not. inp_i(noqp)) goto 8200
          if (.not.rtdb_put(rtdb, 'gw:noqp_beta',mt_int,1,noqp))
     $      call errquit('gw_input: rtdb_put failed', 600, RTDB_ERR)
        elseif (inp_compare(.false.,'vir',test)) then
          if (.not. inp_i(nvqp)) goto 8200
          if (.not.rtdb_put(rtdb, 'gw:nvqp_beta',mt_int,1,nvqp))
     $      call errquit('gw_input: rtdb_put failed', 600, RTDB_ERR)
        else
          goto 8200
        endif
        if (.not. inp_a(test)) goto 10
        goto 602
      endif
c
c     Set first orbital to start counting
c
 700  if (.not.(inp_i(first)))
     $  call errquit('gw_input: "first" must be followed by an integer',
     $                0,INPUT_ERR)
      if (.not.rtdb_put(rtdb,'gw:first',mt_int,1,first))
     $  call errquit('gw_input: rtdb_put failed', 700, RTDB_ERR)
      goto 10 
c
c     Start from NCAP modified eigenvalues
c
 800  if(.not.rtdb_put(rtdb,'gw:ncap',mt_log,1,.true.))
     &  call errquit('gw_input: rtdb_put failed',800, RTDB_ERR)
      goto 10
c
c     Value of imaginary infinitesimal
c
 900  if (.not. inp_f(eta)) 
     &   call errquit('gw_input: eta <float>',0,INPUT_ERR)
      if(.not.rtdb_put(rtdb,'gw:eta',mt_dbl,1,eta))
     &  call errquit('gw_input: rtdb_put failed', 900, RTDB_ERR)
      goto 10
c
c     Compute RPA correlation energy
c
 1000 if(.not.rtdb_put(rtdb,'gw:rpa',mt_log,1,.true.))
     &  call errquit('gw_input: rtdb_put failed',1000, RTDB_ERR)
      goto 10
c
c     Invert matrices by diagonalization 
c
 1100 if(.not.rtdb_put(rtdb,'gw:diag',mt_log,1,.true.))
     &  call errquit('gw_input: rtdb_put failed',1100, RTDB_ERR)
      goto 10
c
c     Start counting from the core upwards
c
 1200 if(.not.rtdb_put(rtdb,'gw:core',mt_log,1,.true.))
     &  call errquit('gw_input: rtdb_put failed',1200, RTDB_ERR)
      goto 10
c
c     normal termination
c
 9999 return
c
c     error handling
c
 8000 call errquit('solver graph || newton <int>',0, INPUT_ERR)
 8200 call errquit('states alpha||beta [occ <int>] [vir <int>]',0,
     &              INPUT_ERR)
c
      end subroutine gw_input    
