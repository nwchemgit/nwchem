      subroutine ccsd_input(rtdb)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 5)
      character*8 dirs(num_dirs)
      character*255 test
      integer maxit, convi, iprt
      data dirs /'maxiter','thresh','        ','print','end'/
c
 10   if (.not. inp_read()) call errquit
     $     ('ccsd_input: inp_read failed', 0)
      if (.not. inp_a(test))
     $     call errquit('ccsd_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('ccsd_input: unknown directive', 0)
c     
      goto (100, 200, 300, 400, 9999) ind
      call errquit('ccsd_input: unimplemented directive', ind)
c     
c     Maximum iterations
c
 100  if (.not. inp_i(maxit)) call errquit
     $     ('ccsd_input: failed reading maxiter', 0)
      if (.not. rtdb_put(rtdb, 'ccsd:maxiter', mt_int, 1, maxit))
     $     call errquit('ccsd-input: rtdb_put failed', 0)
      goto 10
c     
c     thresh <convergence threshold>
c     
 200  if (.not. inp_i(convi)) call errquit
     $     ('ccsd_input: failed reading convergence', 0)
      if (.not. rtdb_put(rtdb, 'ccsd:thresh', mt_int, 1, convi))
     $     call errquit('ccsd-input: rtdb_put failed', 0)
      goto 10
c
c     room for more
c
 300  goto 10
c
c     room for more
c
 400  if (.not. inp_i(iprt)) call errquit
     $     ('ccsd_input: failed reading convergence', 0)
      if (.not. rtdb_put(rtdb, 'ccsd:print', mt_int, 1, iprt))
     $     call errquit('ccsd-input: rtdb_put failed', 0)
      goto 10
c
c     end
c     
 9999 continue
c     
      end

