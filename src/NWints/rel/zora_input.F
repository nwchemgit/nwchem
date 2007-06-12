c     Handles the scalar and spin-orbit ZORA inputs

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
c     <keyword> <scalar||spin-orbit> 
c     
      integer rtdb ! [input] handle to database
c    
      character*255 test

      logical isMolecular
      external isMolecular

c     Preliminaries
      do_zora_first_time = .false.
      do_zora_scalar = .false.
      do_zora_so = .false.
      do_zora_atomic = .true.

c     Read input data beyond zora; store in rtdb.
  10  if (inp_a(test)) then
         if (inp_compare(.false., 'off', test)) then
            do_zora_atomic = .false.

         elseif (inp_compare(.false.,'scalar',test)) then
            do_zora_first_time = .true.
            do_zora_scalar = .true.
            if (isMolecular(test)) do_zora_atomic = .false.

         elseif (inp_compare(.false.,'spin-orbit',test)) then
            do_zora_first_time = .true.
            do_zora_so = .true.
            if (isMolecular(test)) do_zora_atomic = .false.

         else
            call errquit('zora_input: unknown directive',0, UNKNOWN_ERR)

         endif
         goto 10

      endif

czora...debug
c      write(LuOut,*) "do_zora_first_time:",do_zora
c      write(LuOut,*) "do_zora_atomic:",do_zora_atomic
c      write(LuOut,*) "do_zora_scalar:",do_zora_scalar
c      write(LuOut,*) "do_zora_so:",do_zora_so
czora...debug

c
c     Write out a line informing user about the calculation being performed
c
c      if ((ga_nodeid() .eq. 0) .and. do_zora .and. 
c     &     util_print('zora', print_low)) then
c         write(luout,20)
c  20     format(/'ZORA calculations will be performed')
c
c         if (do_zora_scalar) then
c            write(luout,21)
c  21        format(/'Performing scalar ZORA calculations.'/)
c         elseif (do_zora_so) then
c            write(luout,22)
c  22        format(/'Performing spin-orbit ZORA calculations.'/)
c         endif
c      endif
c
c     Put zora parameters in rtdb
c       
      if (.not. rtdb_put(rtdb, 'zora:calc', mt_log, 1, 
     & do_zora_first_time))
     & call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'zora:scalar', mt_log, 1, 
     & do_zora_scalar))
     & call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'zora:so', mt_log, 1, 
     & do_zora_so))
     & call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'zora:atomic', mt_log, 1, 
     & do_zora_atomic))
     & call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)

      return
      end

c     check if we are not using the SAPA method
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
