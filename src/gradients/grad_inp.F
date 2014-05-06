      subroutine gradients_input(rtdb)
C$Id$
c
      implicit none
c
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind,              ! Index of matched directive
     $     checkpoint           ! minutes between checkpoints
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 6)
      character*8 dirs(num_dirs)
      character*255 test
      data dirs /'chkpt','restart','print','noprint','point','end'/
c
 10   if (.not. inp_read()) call errquit
     $     ('gradients_input: inp_read failed', 0, INPUT_ERR)
c
      if (.not. inp_a(test))
     $     call errquit('gradients_input: failed to read keyword', 0,
     &       INPUT_ERR)
      call util_flush(6)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('gradients_input: unknown directive', 0,
     &       INPUT_ERR)

      goto (100, 200, 300, 400, 500, 9999) ind
      call errquit('gradients_input: unimplemented directive', ind,
     &       INPUT_ERR)

 100  if (.not. inp_i(checkpoint)) call errquit
     $     ('scf_input: time between checkpoints missing', 0,
     &       INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'gradients:checkpoint', mt_int, 1, 
     $                   checkpoint))
     $     call errquit('gradients_input: rtdb_put of chkpt failed', 0,
     &       RTDB_ERR)
      goto 10
c     
 200  if (.not. rtdb_put(rtdb,'gradients:restart', MT_LOG, 1, .true.))
     $     call errquit('gradients_input:rtdb_put of restart failed', 0,
     &       RTDB_ERR)
      goto 10
c
c     treatment of point charges now in geometry
c     
 500  call errquit('gradients_input: see geometry bqbq directive',0,
     &       INPUT_ERR)
c
c	print and noprint
c
 300  continue
 400  continue
      call util_print_input(rtdb, 'gradients')
      goto 10
c
c     end
c     
 9999 continue
c     
      end
