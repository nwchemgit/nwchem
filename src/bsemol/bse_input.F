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
      integer, parameter :: num_dirs = 10
      character*11 dirs(num_dirs)
      data dirs / 'method','tda','ncap','ncapr','singlet','triplet', 
     &            'outwin', 'nroots', 'maxroots', 'end' /
c
      integer ind, itemp
      double precision threshold, eta, temp
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
      goto (100, 200, 300, 301, 400, 500, 600, 700, 800, 9999) ind

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
      elseif (inp_compare(.false.,'davidson',test)) then
        if (.not.rtdb_put(rtdb, 'bse:davidson', mt_log, 1, .true.))
     &    call errquit(pname//'rtdb_put failed', 102, RTDB_ERR)
        if (.not.rtdb_put(rtdb, 'bse:analytic', mt_log, 1, .false.))
     &    call errquit(pname//'rtdb_put failed', 101, RTDB_ERR)
      elseif (inp_compare(.false.,'lanczos',test)) then
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
c     use NCAPR DD correction for eigenvalues
c
 301  if (.not.rtdb_put(rtdb, 'bse:ncapr', mt_log, 1, .true.))
     &  call errquit(pname//'rtdb_put failed', 301, RTDB_ERR)
      goto 10
c
c     compute singlet transitions
c
 400  if (.not.rtdb_put(rtdb, 'bse:singlet', mt_log, 1, .true.))
     &  call errquit(pname//'rtdb_put failed', 300, RTDB_ERR)
      if (.not.rtdb_put(rtdb, 'bse:triplet', mt_log, 1, .false.))
     &  call errquit(pname//'rtdb_put failed', 300, RTDB_ERR)
      goto 10
c
c     compute singlet transitions
c
 500  if (.not.rtdb_put(rtdb, 'bse:singlet', mt_log, 1, .false.))
     &  call errquit(pname//'rtdb_put failed', 300, RTDB_ERR)
      if (.not.rtdb_put(rtdb, 'bse:triplet', mt_log, 1, .true.))
     &  call errquit(pname//'rtdb_put failed', 300, RTDB_ERR)
      goto 10
c
c     output energy window
c
 600  if (.not.rtdb_put(rtdb, 'bse:dowindow', mt_log, 1, .true.))
     &  call errquit(pname//'rtdb_put failed', 600, RTDB_ERR)
      if(.not.inp_f(temp))
     &  call errquit(pname//'could not read lower energy',600,INPUT_ERR)
      if(.not.rtdb_put(rtdb,'bse:elower',mt_dbl,1,temp))
     &  call errquit(pname//'rtdb_put failed', 600, RTDB_ERR)
      if(.not.inp_f(temp))
     &  call errquit(pname//'could not read upper energy',600,INPUT_ERR)
      if(.not.rtdb_put(rtdb,'bse:eupper',mt_dbl,1,temp))
     &  call errquit(pname//'rtdb_put failed', 600, RTDB_ERR)
      goto 10
c
c     number of roots
c
 700  if(.not.inp_i(itemp))
     &  call errquit(pname//'could not read no. of roots',700,INPUT_ERR)
      if(.not.rtdb_put(rtdb,'bse:nroots',mt_int,1,itemp))
     &  call errquit(pname//'rtdb_put failed', 700, RTDB_ERR)
      goto 10
c
c     maximum number of trial vectors
c
 800  if(.not.inp_i(itemp))
     &  call errquit(pname//'could not read maxroots',800,INPUT_ERR)
      if(.not.rtdb_put(rtdb,'bse:nmax',mt_int,1,itemp))
     &  call errquit(pname//'rtdb_put failed', 800, RTDB_ERR)
      goto 10
c
c     normal termination
c
 9999 return
c
      end subroutine bse_input
