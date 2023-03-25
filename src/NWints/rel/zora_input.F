c     Handles the spin-free and spin-orbit ZORA inputs

      subroutine zora_input(rtdb)

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
#include "zora.fh"
c
c     
c     Read input for calculations with ZORA
c     
c     All input is of the form:
c     <keyword> <on||off>
c     
      integer rtdb ! [input] handle to database
c    
      character*255 test
c
c     Preliminaries
      do_zora_get_correction = .false.
      do_zora = .false.
c
c     Read input data beyond zora; store in rtdb.
      if (inp_a(test)) then
       if (inp_compare(.false.,'on',test)) then
          do_zora = .true.
       else if (inp_compare(.false.,'off',test)) then
        do_zora = .false.
       else
        call errquit('zora_input: unknown directive',0, UNKNOWN_ERR)
       endif
      else
         call errquit('specify ZORA ON or OFF',0, UNKNOWN_ERR)
      endif
c
c     Put zora parameters in rtdb
      if (.not. rtdb_put(rtdb, 'zora:calc', mt_log, 1,
     & do_zora_get_correction))
     & call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb, 'zora', mt_log, 1,do_zora))
     & call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)
c
      return
      end
c $Id$
