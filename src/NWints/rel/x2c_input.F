c     Handles the X2C input

      subroutine x2c_input(rtdb)

      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rel_consts.fh"
#include "util.fh"
#include "errquit.fh"
#include "x2c.fh"
c
c     
c     Read input for calculations with X2C (exact decoupling)
c     
c     All input is of the form:
c     <keyword> <on||off>
c     
      integer rtdb ! [input] handle to database
c    
      character*255 test
      double precision rtmp
c
c     Preliminaries
      do_x2c = .false.
         
      scale_x2c_hartree = .false.
      scale_x2c_so_hartree = .false.
      x2chscale=1.0d0      ! scaling for Hartree potential
      x2csoscale=1.0d0     ! scaling in the SO operator 
c
c     Read input data beyond x2c; store in rtdb.
  10  if (inp_a(test)) then
       if (inp_compare(.false.,'on',test)) then
        do_x2c = .true.
       else if (inp_compare(.false.,'off',test)) then
         do_x2c = .false.
c
c      process other x2c options
c
c      Hartree scaling?
       else if (inp_compare(.false.,'scalehartree',test)) then
         scale_x2c_hartree = .true.
         rtmp = 1.0d0
         if (.not.inp_f(rtmp)) then          
           write (luout,*) 'rtmp = ',rtmp
           call errquit('x2c: error reading Hartree scaling factor',
     &        1, INPUT_ERR)
         else
           x2chscale = rtmp
         end if

c      Hartree scaling for SO only?
       else if (inp_compare(.false.,'scalesoee',test)) then
         scale_x2c_so_hartree = .true.
         rtmp = 1.0d0
         if (.not.inp_f(rtmp)) then          
           write (luout,*) 'rtmp = ',rtmp
           call errquit('x2c: error reading SOee scaling factor',
     &        1, INPUT_ERR)
         else
           x2csoscale = rtmp
         end if

c
c      keyword not recognized
       else
        call errquit('x2c_input: unknown directive',0, UNKNOWN_ERR)
       endif
       goto 10
      endif
c
c     Put x2c parameters in rtdb

      if (.not. rtdb_put(rtdb, 'x2c:on', mt_log, 1,do_x2c))
     & call errquit('x2c_input: rtdb_put failed', 0, RTDB_ERR)
c
      if (scale_x2c_hartree) then
        write (luout,*) 'X2C: scaling for Hartree term =', x2chscale
      end if
      if (scale_x2c_so_hartree) then
        write (luout,*) 'X2C: scaling for SO Hartree term =',
     &     x2csoscale
      end if
c
      return
      end
