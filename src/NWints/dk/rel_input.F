      subroutine rel_input(rtdb)
*
* $Id: rel_input.F,v 1.6 2000-03-16 22:30:53 mg201 Exp $
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
c     Read input for the relativistic block
c     
      integer rtdb ! [input] handle to database
c     
      character*255 test
      double precision clight
      logical dk_run, dmd_run, status
c
c     Set the two relativistic approaches to .false.
c
      dk_run = .false.
      if (.not. rtdb_put(rtdb, 'doug_kroll', mt_log, 1, dk_run))
     $   call errquit('rel_input: rtdb_put failed', 0)
      dmd_run = .false.
      if (.not. rtdb_put(rtdb, 'dyall_mod_dir', mt_log, 1, dmd_run))
     $   call errquit('rel_input: rtdb_put failed', 0)
c
c     Set some relativistic variables we could need later on
c
   1  if (.not. inp_read()) call errquit('rel_input: inp_read failed',0)
      if (.not. inp_a(test))
     $     call errquit('rel_input: failed to read keyword', 0)
c
      if (inp_compare(.false.,'clight', test)) then
c
         if (.not. inp_f(clight)) then
            call errquit('rel_input: failed reading clight', 0)
         else 
           if(.not.rtdb_put(rtdb,'relativistic:clight',mt_dbl,1,clight))
     $        call errquit('rel_input: rtdb_put failed', 0)
         endif
         write(luout,10) clight
  10        format(/,'The speed of light used in the relativistic ',
     $               'calculations is ',F16.10, 'au.'/)
c
      elseif (inp_compare(.false.,'douglas-kroll', test)) then
c
         call dk_input(rtdb)
c
      elseif (inp_compare(.false.,'dyall-mod-dirac', test)) then
c
         call dmd_input(rtdb)
c
      else if (inp_compare(.false.,'end', test)) then
c
c        end will be catched in while loop structure
c
      else
         call errquit('rel_input: unrecognized directive', 0)
      endif
c
c     End of while loop over the input
c
      if (.not. inp_compare(.false.,'end', test)) goto 1
c
c     User could have set both methods by accident. Check if this is the case 
c     stop if both methods are set to true.
c
      status = rtdb_get(rtdb, 'doug_kroll', mt_log, 1, dk_run)
      status = rtdb_get(rtdb, 'dyall_mod_dir', mt_log, 1, dmd_run)
      if (dk_run .and. dmd_run) then
         if (ga_nodeid() .eq. 0) then
         write(luout,20)
  20     format(/'Input error relativistic:'/,
     $   'Both Douglas-Kroll and Dyall-Mod-Dirac have been specified !')
         endif
         call errquit('rel_input: input error',0)
      endif
c
      return
c       
      end
