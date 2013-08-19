c $Id$
      subroutine et_input(rtdb)
      implicit none
#include "errquit.fh"
c
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      integer rtdb, method_2e, fmo
      character*255 test
      double precision tol2e
c
 100  if (.not. inp_read()) call errquit('et_input: inp_read failed',0,
     &       INPUT_ERR)
c
       if (.not. inp_a(test))
     $     call errquit('scf_input: failed to read keyword', 0,
     &       INPUT_ERR)
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
     $        ('et_input: tol2e <thresh>',0, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'et:tol2e', MT_DBL, 1, tol2e))
     $        call errquit('et_input: rtdb_put of tol2e failed',0,
     &       RTDB_ERR)
c 
      goto 100
c
c     choice of 2e method is read
c     ---------------------------
      else if (inp_compare(.false.,'fock', test)) then
        method_2e = 1
         if (.not. rtdb_put(rtdb,'et:method_2e', MT_INT, 1, method_2e))
     $        call errquit('et_input: rtdb_put of method_2e failed',0,
     &       RTDB_ERR)
c
      goto 100
c
      else if (inp_compare(.false.,'nofock', test)) then
        method_2e = 2
c
         if (.not. rtdb_put(rtdb,'et:method_2e', MT_INT, 1, method_2e))
     $        call errquit('et_input: rtdb_put of method_2e failed',0,
     &       RTDB_ERR)
c
      goto 100
c
c     choice of FMO method
c     ---------------------------
      else if (inp_compare(.false.,'fmo', test)) then
        fmo = 1
         if (.not. rtdb_put(rtdb,'et:fmo', MT_INT, 1, fmo))
     $        call errquit('et_input: rtdb_put of fmo failed',0,
     &       RTDB_ERR)
c
      goto 100
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
