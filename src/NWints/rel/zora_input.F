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
c     Read input for calculations with ZORA
c     
c     All input is of the form:
c     <keyword> <on||off>          or
c     <keyword> <sf||so> 
c     
      integer rtdb ! [input] handle to database
c    
      character*255 test

      logical isMolecular
      external isMolecular

c     Preliminaries
      do_zora_get_correction = .false.
      do_zora_sf = .false.
      do_zora_so = .false.
      do_zora = .false.

c     Read input data beyond zora; store in rtdb.
  10  if (inp_a(test)) then
       if (inp_compare(.false.,'sf',test)) then
            do_zora_get_correction = .true.
            do_zora_sf = .true.
       elseif (inp_compare(.false.,'so',test)) then
            do_zora_get_correction = .true.
            do_zora_so = .true.
       else
            call errquit('zora_input: unknown directive',0, UNKNOWN_ERR)
       endif
       goto 10
      endif
c
c     Write out a line informing user about the calculation being performed
c
      if (do_zora_sf) then
        write(luout,21)
  21    format(/'Performing spin-free ZORA calculations.'/)
      elseif (do_zora_so) then
        write(luout,22)
  22    format(/'Performing spin-orbit ZORA calculations.'/)
      endif
c
c     Put zora parameters in rtdb
c       
      if (.not. rtdb_put(rtdb, 'zora:calc', mt_log, 1, 
     & do_zora_get_correction))
     & call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'zora:sf', mt_log, 1,do_zora_sf))
     & call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'zora:so', mt_log, 1,do_zora_so))
     & call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)

      return
      end

c
c     check if we are not using the SAPA method
c
      logical function isMolecular(test)

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

      character*255 test

c     Preliminaries
      isMolecular = .false.

c     Is this is a self-consistent calculation ?
      if (inp_a(test)) then

        if (inp_compare(.false., 'molecular', test)) then
          isMolecular = .true.
        else if (inp_compare(.false., 'atomic', test)) then
          isMolecular = .false.
        end if

      end if

      return
      end 
