c $Id: et_input.F,v 1.4 2002-02-19 21:38:23 d3k958 Exp $
      subroutine et_input(rtdb)
      implicit none
c
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      integer rtdb
      character*255 test
c
 100  if (.not. inp_read()) call errquit('et_input: inp_read failed',0)
c
       if (.not. inp_a(test))
     $     call errquit('scf_input: failed to read keyword', 0)
c     
c     print level is read
c     -------------------
      if (inp_compare(.false.,'print', test)) then
         call util_print_input(rtdb, "et")
c
      goto 100
c
c     MO vector filenames are read
      else if (inp_compare(.false.,'vectors', test)) then
         call et_vecs_input(rtdb,'et')
c 
         goto 100
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
