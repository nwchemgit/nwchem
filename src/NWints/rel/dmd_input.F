      subroutine dmd_input(rtdb)
*
* $Id: dmd_input.F,v 1.1 1999-12-16 23:44:05 bert Exp $
*
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "rel_consts.fh"
c     
c     Read input for calculations with Dyall-Modified-Dirac or Dyall-Mod-Dirac
c     
c     All input is of the form:
c     <keyword>                          or   
c     <keyword> <on||off>                or
c     <keyword> <sf>  !future option     or
c     <keyword> <afw> !future option     or a combination of the last three
c     
      integer rtdb ! [input] handle to database
c    
      character*255 test
      logical dmd_run
c
c     Read input data beyond Dyall-Modified-Dirac or Dyall-Mod-Dirac; store in rtdb.
c     
  10  if (inp_a(test)) then
         if (inp_compare(.false., 'off', test)) then
            dmd_run = .false.
         elseif (inp_compare(.false., 'on', test)) then
            dmd_run  = .true.
         else
            call errquit('dk_input: unknown directive',0)
         endif
         goto 10
      endif
c
c     If we are going to do Dyall-Modified-Dirac then print out a line informing user,
c
      if ((ga_nodeid() .eq. 0) .and . dmd_run .and. 
     $     util_print('dyall-mod-dirac', print_low)) then
         write(luout,20)
  20     format(/'Explain that this is Dyall-Modified-Dirac !')
      endif
c
c     Put parameters in rtdb
c       
      if (.not. rtdb_put(rtdb, 'dyall_mod_dir', mt_log, 1, dmd_run))
     $   call errquit('dk_input: rtdb_put failed', 0)
c
      return
c       
      end
