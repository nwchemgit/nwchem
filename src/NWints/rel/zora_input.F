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
      integer zora_type
      logical zora_run

      logical isMolecular
      external isMolecular

c     Preliminaries
      zora_run = .true.
      do_zora = .true.
      do_scalar_zora = .false.
      do_so_zora = .false.
      zora_atomic = 1  ! default: use sum of atomic potential approximation (SAPA)

c     Read input data beyond zora; store in rtdb.
  10  if (inp_a(test)) then
         if (inp_compare(.false., 'off', test)) then
            zora_type = 0
            zora_run = .false.
            do_zora = .false.
            return

         elseif (inp_compare(.false.,'on',test)) then
            zora_type = 1   ! perform scalar zora by default
            do_scalar_zora = .true.
            if (isMolecular(test)) zora_atomic = 0

         elseif (inp_compare(.false.,'scalar',test)) then
            zora_type = 1
            do_scalar_zora = .true.
            if (isMolecular(test)) zora_atomic = 0

         elseif (inp_compare(.false.,'spin-orbit',test)) then
            zora_type = 2
            do_so_zora = .true.
            if (isMolecular(test)) zora_atomic = 0
            call errquit('zora_input:
     & Spin-orbit ZORA not implemented yet.',0, UNKNOWN_ERR)

         else
            call errquit('zora_input: unknown directive',0, UNKNOWN_ERR)

         endif
         goto 10

      endif

cdebug
c      write(LuOut,*) "do_zora:",do_zora
c      write(LuOut,*) "zora_atomic:",zora_atomic
cdebug

c
c     Write out a line informing user about the calculation being performed
c
      if ((ga_nodeid() .eq. 0) .and. zora_run .and. 
     $     util_print('zora', print_low)) then
         write(luout,20)
  20     format(/'ZORA calculations will be performed')

         if (zora_type .eq. 1) then
            write(luout,21)
  21        format(/'Performing scalar ZORA calculations.'/)
         elseif (zora_type .eq. 2) then
            write(luout,22)
  22        format(/'Performing spin-orbit ZORA calculations.'/)

         endif

      endif
c
c     Put parameters in rtdb
c       
      if (.not. rtdb_put(rtdb, 'zora', mt_log, 1, zora_run))
     $   call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'zora:type', mt_int, 1, zora_type))
     $   call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb, 'zora:atomic', mt_int, 1, zora_atomic))
     $   call errquit('zora_input: rtdb_put failed', 0, RTDB_ERR)

      return
      end


c     check if we are not using the APA method
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
        else
          isMolecular = .false.
        end if
      end if

      return
      end 
