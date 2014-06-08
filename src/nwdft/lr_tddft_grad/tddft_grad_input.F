c
c     tddft gradients input handler
c
      subroutine tddft_grad_input(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "geom.fh"
#include "stdio.fh"
#include "case.fh"
c
      integer rtdb
c
      integer iroot                 ! specific root
      double precision solve_thresh ! convergence threshold for solver
      integer solve_maxit           ! maximum number of iterations in solver
      integer solve_maxsub          ! maximum number expansion vectors per root
c
      character*255 test
      character*32 pname
c
      pname="tddft_grad_input: "
c
c     --------
c     Defaults
c     --------
c
c     TDDFT gradient settings
c
      iroot = 1
      if (.not.rtdb_put(rtdb,'tddft_grad:isinglet_roots',mt_int,1,
     & iroot)) call errquit(pname//'failed to write iroot',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'tddft_grad:itriplet_roots',mt_int,1,
     & iroot)) call errquit(pname//'failed to write iroot',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb,'tddft_grad:iroots',mt_int,1,iroot))
     &  call errquit(pname//'failed to write iroot',0,RTDB_ERR)
c
      solve_thresh = 1d-04
      if (.not. rtdb_put(rtdb,'tddft_grad:solve_thresh',mt_dbl,
     &           1,solve_thresh))
     &  call errquit(pname//'failed to write solve_thresh',0,RTDB_ERR)
c
c    Default consistent with what is present 
c    in the TDDFT gradient routine (tddft_grad_init.F)
c
      solve_maxsub = 40
      if (.not. rtdb_put(rtdb,'tddft_grad:solve_maxsub',mt_int,
     &           1,solve_maxsub))
     &  call errquit(pname//'failed to write solve_maxsub',0,RTDB_ERR)
c 
      solve_maxit = 250
      if (.not. rtdb_put(rtdb,'tddft_grad:solve_maxit',mt_int,
     &           1,solve_maxit))
     &  call errquit(pname//'failed to write solve_maxit',0,RTDB_ERR)
c
c     ----------
c     Read input
c     ----------
c
 10   if (.not. inp_read())
     1  call errquit(pname//'failed reading input',0,INPUT_ERR)
      if (.not. inp_a(test))
     1  call errquit(pname//'failed reading keyword',0,INPUT_ERR)
c
c     Which root ?
c
      if (inp_compare(.false.,test,'root')) then
       if (.not.inp_i(iroot)) iroot = 1
      if (.not. rtdb_put(rtdb,'tddft_grad:isinglet_roots',mt_int,1,
     & iroot)) call errquit(pname//'failed to write iroot',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb,'tddft_grad:itriplet_roots',mt_int,1,
     & iroot)) call errquit(pname//'failed to write iroot',0,RTDB_ERR)
       if (.not. rtdb_put(rtdb,'tddft_grad:iroots',mt_int,1,iroot))
     &  call errquit(pname//'failed to write iroot',0,RTDB_ERR)
c
c     Solver threshold
c
      else if (inp_compare(.false.,test,'solve_thresh')) then
       if (.not.inp_f(solve_thresh)) solve_thresh = 1.d-04
       if (.not. rtdb_put(rtdb,'tddft_grad:solve_thresh',mt_dbl,
     &           1,solve_thresh))
     &  call errquit(pname//'failed to write solve_thresh',0,RTDB_ERR)
c
c     Solver subspace
c
      else if (inp_compare(.false.,test,'solve_maxsub')) then
       if (.not.inp_i(solve_maxsub)) solve_maxsub = 10
       if (.not. rtdb_put(rtdb,'tddft_grad:solve_maxsub',mt_int,
     &           1,solve_maxsub))
     &  call errquit(pname//'failed to write solve_maxsub',0,RTDB_ERR)
c 
c     Solver iterations
c
      else if (inp_compare(.false.,test,'solve_maxit')) then
       if (.not.inp_i(solve_maxit)) solve_maxit = 250
       if (.not. rtdb_put(rtdb,'tddft_grad:solve_maxit',mt_int,
     &           1,solve_maxit))
     &  call errquit(pname//'failed to write solve_maxit',0,RTDB_ERR)
c
c     Print
c
      else if (inp_compare(.false.,test,'print')) then
        call util_print_input(rtdb,'tddft_grad')
c
c     end
c
      else if (inp_compare(.false.,test,'end')) then
        goto 20
c
      else
        call errquit(pname//'unknown directive',0, INPUT_ERR)
      endif
      goto 10
c
 20   return
      end
c $Id$
