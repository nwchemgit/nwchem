      subroutine rel_input(rtdb)
*
* $Id: rel_input.F,v 1.1 1999-12-15 22:46:48 bert Exp $
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
c     Read input for the relativistic block
c     
      integer rtdb ! [input] handle to database
c     
      character*255 test
      double precision clight
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
           if (.not. rtdb_put(rtdb,'doug_kroll:clight',mt_dbl,1,clight))
     $        call errquit('rel_input: rtdb_put failed', 0)
         endif
c
      elseif (inp_compare(.false.,'douglas-kroll', test)) then
c
         call dk_input(rtdb)
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
      return
c       
      end
