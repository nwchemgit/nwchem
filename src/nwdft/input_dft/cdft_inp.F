      subroutine cdft_input(rtdb, nconstr, dtype, iatms1, iatms2, 
     $                      fatms1, fatms2, constr)
c
c     $Id$
c
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
#include "cdft.fh"
      integer rtdb ! [input] handle to database
c
      integer ind               ! Index of matched directive
      integer num_dirs, num_c, num_p         ! No. of known directives
      parameter (num_dirs = 2, num_c=2, num_p=3) 
      character*12 dirs(num_dirs), dirs_c(num_c), dirs_p(num_p)
      character*255 test
      integer nconstr
      integer dtype(max_constr), iatms1(max_constr), fatms1(max_constr)
      integer iatms2(max_constr), fatms2(max_constr)
      double precision constr(max_constr)
      double precision dl_conv  
      logical grid_make_sa
      integer ipop
      data dirs /'convergence', 'pop'/
      data dirs_c /'charge', 'spin'/
      data dirs_p /'becke', 'mulliken', 'lowdin'/
c
c
      if (.not.rtdb_get(rtdb,'dft:grid_make_sa',mt_log,1,grid_make_sa))
     &   grid_make_sa = .false.

      if (.not.(inp_i(iatms1(nconstr)).and.inp_i(fatms1(nconstr))))
     &   call errquit('cdft_input: atom not selected', 2700, RTDB_ERR)
      if (.not.(inp_i(iatms2(nconstr)).and.inp_i(fatms2(nconstr)))) then
         iatms2(nconstr) = 0
         fatms2(nconstr) = 0
      endif
      if (.not. inp_a(test)) call errquit
     $         ('cdft_input: specify charge or spin', 2700, RTDB_ERR)
      if (.not. inp_match(num_c, .false., test, dirs_c, ind))
     $     call errquit('cdft_input: unknown constraint', 0, INPUT_ERR)

      goto (1, 2, 19) ind
 1    dtype(nconstr) = 1
      goto 19
 2    dtype(nconstr) = 2
      goto 19
 19   continue
      if (.not. inp_f(constr(nconstr)))
     &  call errquit('cdft: constraint value not defined',2700,RTDB_ERR)

 10   if (.not. inp_a(test)) goto 1999
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('cdft_input: unknown directive', 0,
     &       INPUT_ERR)
c     
      goto (100, 200, 1999) ind 
      call errquit('cdft_input: unimplemented directive', ind,
     &       INPUT_ERR)

 100  if (.not. inp_f(dl_conv))
     &  dl_conv = 1d-6
      if (.not. rtdb_put(rtdb,'dft:dl_conv', mt_dbl, 1, dl_conv))
     &   call errquit('cdft_input: rtdb_put failed', 2500, RTDB_ERR)
      goto 10

 200  if (.not. inp_a(test)) then
       if (.not. rtdb_get(rtdb,'dft:ipop', mt_int, 1, ipop))
     &   ipop = 3
        goto 299
      endif

      if (.not. inp_match(num_p, .false., test, dirs_p, ind))
     $     call errquit('cdft: unknown population', 0, INPUT_ERR)
c
      goto (201, 202, 203, 299) ind

 201  ipop = 1
      grid_make_sa = .true.
      goto 299
 202  ipop = 2
      goto 299
 203  ipop = 3
      goto 299
 299  continue
      if (.not. rtdb_put(rtdb,'dft:ipop', mt_int, 1, ipop))
     &   call errquit('cdft_input: rtdb_put failed', 2800, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'dft:grid_make_sa',mt_log,1,grid_make_sa))
     &   call errquit('cdft_input: rtdb_put failed', 2900, RTDB_ERR)
      goto 10

 1999 continue
      return 
      end
