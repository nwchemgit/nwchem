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
c     <keyword> <scalar||spin-orbit> or a combination of the last two
c     
      integer rtdb ! [input] handle to database
c    
      character*255 test
      integer zora_type
      logical zora_run

      zora_run = .true.
c
c     Read input data beyond zora; store in rtdb.
c     
  10  if (inp_a(test)) then
         if (inp_compare(.false., 'off', test)) then
            zora_type = 0
            zora_run = .false.
         elseif (inp_compare(.false.,'on',test)) then
            zora_type = 1   ! perform scalar zora by default
         elseif (inp_compare(.false.,'scalar',test)) then
            zora_type = 1
         elseif (inp_compare(.false.,'spin-orbit',test)) then
            zora_type = 2
            call errquit('zora_input:
     & Spin-orbit ZORA not implemented yet.',0, UNKNOWN_ERR)
         else
            call errquit('zora_input: unknown directive',0, UNKNOWN_ERR)
         endif
         goto 10
      endif
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

      return
      end
