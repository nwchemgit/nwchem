c $Id: et_input.F,v 1.5 2002-03-15 15:55:37 d3k958 Exp $
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
      double precision toldens, tol2e
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
c     ----------------------------
      else if (inp_compare(.false.,'vectors', test)) then
         call et_vecs_input(rtdb,'et')
c 
         goto 100
c
c     screening tolerances are read
c     -----------------------------
      else if (inp_compare(.false.,'tol2e', test)) then
         if (.not. inp_f(tol2e)) call errquit
     $        ('et_input: tol2e <thresh>',0)
         if (.not. rtdb_put(rtdb,'et:tol2e', MT_DBL, 1, tol2e))
     $        call errquit('et_input: rtdb_put of tol2e failed',0)
c 
         goto 100
c
      else if (inp_compare(.false.,'toldens', test)) then
         if (.not. inp_f(toldens)) call errquit
     $        ('et_input: toldens <thresh>',0)
         if (.not. rtdb_put(rtdb,'et:toldens', MT_DBL, 1, toldens))
     $        call errquit('et_input: rtdb_put of tol2e failed',0)
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
