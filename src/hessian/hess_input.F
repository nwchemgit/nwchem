      subroutine hess_input(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "nwc_const.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer rtdb
      double precision conv
      character*255 field
c
c     ----- read values from input -----
c
 10   if(.not.inp_read())
     $     call errquit('hess_input: unexpected eof', 0, INPUT_ERR)
      if (.not. inp_a(field))
     $     call errquit('hess_input: failed to read field', 0,
     &       INPUT_ERR)
c
      if(inp_compare(.false.,'end',field)) then
         return
c
      else if(inp_compare(.false.,'print',field) .or.
     $        inp_compare(.false.,'noprint',field)) then
         call util_print_input(rtdb,'hess')
c
      else if(inp_compare(.false.,'profile',field)) then
         if (.not. rtdb_put(rtdb,'hess:pstat', MT_LOG, 1, .true.))
     $      call errquit('hess_input:rtdb_put of pstat failed', 0,
     &       INPUT_ERR)
c
      else if (inp_compare(.false.,'thresh',field)) then
         if (.not. inp_f(conv)) call errquit
     $        ('hess_input: failed reading convergence', 0,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'hessian:thresh', mt_dbl, 1, conv))
     $        call errquit('scf-input: rtdb_put failed', 0, RTDB_ERR)
      else
         call errquit('hess_input: unknown directive',0, INPUT_ERR)
      endif
c
      goto 10                   ! While loop
c
      end
