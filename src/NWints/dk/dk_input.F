      subroutine dk_input(rtdb)
*
* $Id: dk_input.F,v 1.5 1999-12-15 22:45:51 bert Exp $
*
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
c     
c     Read input for calculations with Douglas-Kroll
c     
c     All input is of the form:
c     <keyword>                    or   
c     <keyword> <on||off>          or
c     <keyword> <fpp||dkh||dkfull> or a combination of the last two
c     
      integer rtdb ! [input] handle to database
      
      character*255 test
      logical dk_run, status
      integer dk_type
c
c     Default settings
c
      dk_run    = .true.
      dk_type   = 2
c
c     Check for a previous DK calculation and if the type was set there.
c
      status = rtdb_get(rtdb, 'doug_kroll:type', MT_INT, 1,dk_type)
c
c     Read input data beyond douglas-kroll; store in rtdb
c     
  10  if (inp_a(test)) then
         if (inp_compare(.false., 'off', test)) then
            dk_run  = .false.
         elseif (inp_compare(.false., 'on', test)) then
            dk_run  = .true.
         elseif (inp_compare(.false.,'fpp',test)) then
            dk_type = 1
         elseif (inp_compare(.false.,'dkh',test)) then
            dk_type = 2 
         elseif (inp_compare(.false.,'dkfull',test)) then
            dk_type = 3
         else
            call errquit('dk_input: unknown directive',0)
         endif
         goto 10
      endif
c
c     If we are going to do Douglas-Kroll, write out a line informing user
c
      if ((ga_nodeid() .eq. 0) .and .dk_run .and. 
     $     util_print('douglas-kroll', print_low)) then
         write(luout,20)
  20     format(/12x, 'Douglas-Kroll integrals are switched on !'/)
         if (dk_type .eq. 1) then
            write(luout,21)
  21        format(/12x, 'Integrals for the no-pair equation with ',
     $             'free-particle projection operators'/
     $             12x,'will be computed.'//)
         elseif (dk_type .eq. 2) then
            write(luout,22)
  22        format(/12x, 'Integrals for the no-pair equation with ',
     $             'external-field projection operators'/
     $             12x,'will be computed (Douglas-Kroll-Hess).'//)
         else  
            write(luout,23)
  23        format(/12x, 'Integrals for the no-pair equation with ',
     $             'external-field projection operators'/12x,
     $             'including cross-product terms will be computed.'//)
         endif
         call util_flush(luout)
      endif
c
c     Put parameters in rtdb
c       
      if (.not. rtdb_put(rtdb, 'doug_kroll', mt_log, 1, dk_run))
     $   call errquit('dk_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'doug_kroll:type', mt_int, 1, dk_type))
     $   call errquit('dk_input: rtdb_put failed', 0)
        return
c       
        end
