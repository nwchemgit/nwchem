
      subroutine string_input(rtdb)

      implicit none
! External stuff
#include "errquit.fh"
#include "inp.fh"
#include "nwc_const.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
c   passed in variable
c
      integer rtdb
c
c   local variable
c
      character*255 field
      real*8 tol, temp, stepsize
      integer nbeads, maxit, interpol
      logical frze
c
c     ----- read values from input -----
c
 10   if(.not.inp_read())
     $     call errquit('string_input: unexpected eof',911, INPUT_ERR)
      if (.not. inp_a(field))
     $     call errquit('string_input: failed to read field',911,
     &       INPUT_ERR)
c
      if(inp_compare(.false.,'end',field)) then
         return
c
      else if(inp_compare(.false.,'tol',field)) then
         if(.not. inp_f(tol))
     $        call errquit('string_input: tol',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'string:tol',mt_dbl,1,tol))
     $        call errquit('string_input: rtdb put failed',0,
     &       RTDB_ERR)
c
      else if(inp_compare(.false.,'temperature',field)) then
         if(.not. inp_f(temp))
     $        call errquit('string_input: temperature',911,
     &       INPUT_ERR)
         if (temp .ne. 0)
     $        call errquit('Only zero temperature implemented',0,
     &       INPUT_ERR)
c
      else if(inp_compare(.false.,'stepsize',field)) then
         if(.not. inp_f(stepsize))
     $        call errquit('string_input: stepsize',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'string:stepsize',mt_dbl,1,stepsize))
     $        call errquit('string_input: rtdb put failed',0,
     &       RTDB_ERR)
c
      else if (inp_compare(.false.,'nbeads', field)) then
         if(.not. inp_i(nbeads))
     $        call errquit('string_input: nbeads',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'string:nbeads',mt_int,1,
     $        nbeads))
     $        call errquit('string_input: rtdb put failed',0, RTDB_ERR)
         if (nbeads .lt. 3)
     $        call errquit('string_input: nbeads less than 3',911,
     $        INPUT_ERR)
c
      else if (inp_compare(.false.,'maxiter', field)) then
         if(.not. inp_i(maxit))
     $        call errquit('string_input: maxiter',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'string:maxit',mt_int,1,maxit))
     $        call errquit('string_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'interpol', field)) then
         if(.not. inp_i(interpol))
     $        call errquit('string_input: interpol',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'string:interpol',mt_int,1,
     $        interpol))
     $        call errquit('string_input: rtdb put failed',0, RTDB_ERR)

      else if (inp_compare(.false.,'freeze1', field)) then
         if (.not. inp_l(frze)) frze = .false.
         if (frze) then
           if (.not. rtdb_put(rtdb,'string:freeze1',mt_log,1,frze))
     $        call errquit('string_input: rtdb put failed',0, RTDB_ERR)
         end if

      else if (inp_compare(.false.,'freezen', field)) then
         if (.not. inp_l(frze)) frze = .false.
         if (frze) then
           if (.not. rtdb_put(rtdb,'string:freezen',mt_log,1,frze))
     $        call errquit('string_input: rtdb put failed',0, RTDB_ERR)
         end if
c
      else
         call errquit('string_input: unknown directive',0, INPUT_ERR)
      endif
c
      goto 10                   ! While loop
c
      end

c $Id$
