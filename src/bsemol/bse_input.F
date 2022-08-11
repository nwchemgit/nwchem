      subroutine bse_input(rtdb)
c
c     Input parser for molecular BSE module
c
c     Daniel Mejia-Rodriguez
c     November 2022
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
      character(*), parameter :: pname = "bse_input: "
c
      integer, parameter :: num_dirs = 4
      character*11 dirs(num_dirs)
      data dirs / 'method', 'tda', 'ncap', 'end' /
c
      integer ind, qpiter,eviter,noqp(2),nvqp(2), ngl, first
      double precision threshold, eta
      double precision, parameter :: ha2ev = 27.211386245988d0
c
c     Set default values and put them in RTDB
c
      call bse_defaults(rtdb)
c
c     Parse the input file
c
 10   if (.not. inp_read()) then
        call errquit(pname//'failed to read input',0,INPUT_ERR)
      endif
c
      if (.not. inp_a(test)) then
        call errquit(pname//'failed to read keyword',0,INPUT_ERR)
      endif
c
      if (.not.inp_match(num_dirs, .false., test, dirs, ind)) then
        call errquit(pname//'unknown directive', 0, INPUT_ERR)
      end if
c
      goto (100, 200, 300, 9999) ind

      call errquit(pname//'unknown directive',ind, INPUT_ERR)
c
c     -----------------------------------------------------------------
c
c     method
c
 100  if (.not. inp_a(test)) 
     &  call errquit(pname//'could not read method', 100, INPUT_ERR)

      if (inp_compare(.false.,'analytic',test)) then
        if (.not.rtdb_put(rtdb, 'bse:analytic', mt_log, 1, .true.))
     &    call errquit(pname//'rtdb_put failed', 101, RTDB_ERR)
      elseif (inp_compare(.false.,'bse:davidson',test)) then
        if (.not.rtdb_put(rtdb, 'bse:davidson', mt_log, 1, .true.))
     &    call errquit(pname//'rtdb_put failed', 102, RTDB_ERR)
        if (.not.rtdb_put(rtdb, 'bse:analytic', mt_log, 1, .false.))
     &    call errquit(pname//'rtdb_put failed', 101, RTDB_ERR)
      elseif (inp_compare(.false.,'bse:lanczos',test)) then
        if (.not.rtdb_put(rtdb, 'bse:lanczos', mt_log, 1, .true.))
     &    call errquit(pname//'rtdb_put failed', 103, RTDB_ERR)
        if (.not.rtdb_put(rtdb, 'bse:analytic', mt_log, 1, .false.))
     &    call errquit(pname//'rtdb_put failed', 101, RTDB_ERR)
      else
        call errquit(pname//'unknown BSE method',104,INPUT_ERR)
      endif
      goto 10
c
c     use Tamm-Dancoff approximation
c
 200  if (.not.rtdb_put(rtdb, 'bse:tda', mt_log, 1, .true.))
     &  call errquit(pname//'rtdb_put failed', 200, RTDB_ERR)
      goto 10
c
c     use NCAP DD correction for eigenvalues
c
 300  if (.not.rtdb_put(rtdb, 'bse:ncap', mt_log, 1, .true.))
     &  call errquit(pname//'rtdb_put failed', 300, RTDB_ERR)
      goto 10
c
c     normal termination
c
 9999 return
c
      end subroutine bse_input
