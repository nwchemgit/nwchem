      subroutine et_input(rtdb)
      implicit none
c
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      integer rtdb
      logical status
      character*255 test
c
c
c     
 1     if (.not. inp_read()) call errquit('et_input: inp_read failed',0)
c
      if (.not. inp_a(test))
     $     call errquit('scf_input: failed to read keyword', 0)
c     
c     print level is read
c     -------------------
      if (inp_compare(.false.,'print', test)) then
         call util_print_input(rtdb, "et")
         goto 1
c
c
c     end of input block
c     -------------------
      else if (.not. inp_compare(.false.,'end', test)) then
       continue
      endif
c     
c
      return
      end
