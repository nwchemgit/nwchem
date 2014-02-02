c
c     gradients block input handler
c
      subroutine gradients_input(rtdb)
c
C$Id$
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
      integer checkpoint   ! minutes between checkpoints
      logical restart
c
      integer iroot
c
      double precision solve_thresh
      integer solve_maxsub
      integer solve_maxit
c
      character*255 test
      character*32 pname
c
      pname="gradients_input: "
c
c     --------
c     Defaults
c     --------
c
      checkpoint = 2
      if (.not. rtdb_put(rtdb,'gradients:checkpoint',mt_int,1,
     &          checkpoint))
     &  call errquit(pname//'failed to write checkpoint',0,RTDB_ERR)
c
      restart = .false.
       if (.not. rtdb_put(rtdb,'gradients:restart',mt_log,1,.true.))
     &  call errquit(pname//'failed to write restart',0,RTDB_ERR)
c
c     For excited-state gradient calculations
c
      iroot = 1
      if (.not. rtdb_put(rtdb,'gradients:isinglet_roots',mt_int,1,
     & iroot)) call errquit(pname//'failed to write iroot',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb,'gradients:itriplet_roots',mt_int,1,
     & iroot)) call errquit(pname//'failed to write iroot',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb,'gradients:iroots',mt_int,1,iroot))
     &  call errquit(pname//'failed to write iroot',0,RTDB_ERR)
c
      solve_thresh = 1d-04
      if (.not. rtdb_put(rtdb,'gradients:solve_thresh',mt_dbl,
     &           1,solve_thresh))
     &  call errquit(pname//'failed to write solve_thresh',0,RTDB_ERR)
c
c    Default consistent with what is present 
c    in the TDDFT gradient routine (src/nwdft/lr_tddft_grad/tddft_grad_init.F).
      solve_maxsub = 40
      if (.not. rtdb_put(rtdb,'gradients:solve_maxsub',mt_int,
     &           1,solve_maxsub))
     &  call errquit(pname//'failed to write solve_maxsub',0,RTDB_ERR)
c 
      solve_maxit = 250
      if (.not. rtdb_put(rtdb,'gradients:solve_maxit',mt_int,
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
c     Minutes between checkpoints otal number of nuclear time steps: nstep_nucl
      if (inp_compare(.false.,test,'checkpoint')) then
        if (.not.inp_i(checkpoint)) then
          write(luOut,*) 'checkpoint: value not found; ',
     1      'default value of 2 minutes will be used'
          checkpoint = 2
        endif
        if (.not. rtdb_put(rtdb,'gradients:checkpoint',mt_int,1,
     &          checkpoint))
     &  call errquit(pname//'failed to write checkpoint',0,RTDB_ERR)
c
c     Restart
      else if (inp_compare(.false.,test,'restart')) then
       if (.not. rtdb_put(rtdb,'gradients:restart',mt_log,1,.true.))
     &  call errquit(pname//'failed to write restart',0,RTDB_ERR)
c
c     Point
      else if (inp_compare(.false.,test,'point')) then
       call errquit(pname//'see geometry bqbq directive',0,INPUT_ERR)
c
c     For excited-state gradient calculations
c
c     Which root ?
      else if (inp_compare(.false.,test,'root')) then
       if (.not.inp_i(iroot)) iroot = 1
      if (.not. rtdb_put(rtdb,'gradients:isinglet_roots',mt_int,1,
     & iroot)) call errquit(pname//'failed to write iroot',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb,'gradients:itriplet_roots',mt_int,1,
     & iroot)) call errquit(pname//'failed to write iroot',0,RTDB_ERR)
       if (.not. rtdb_put(rtdb,'gradients:iroots',mt_int,1,iroot))
     &  call errquit(pname//'failed to write iroot',0,RTDB_ERR)
c
c     Linear system solver settings
      else if (inp_compare(.false.,test,'solve_thresh')) then
       if (.not.inp_f(solve_thresh)) solve_thresh = 1.d-04
       if (.not. rtdb_put(rtdb,'gradients:solve_thresh',mt_dbl,
     &           1,solve_thresh))
     &  call errquit(pname//'failed to write solve_thresh',0,RTDB_ERR)
c
      else if (inp_compare(.false.,test,'solve_maxsub')) then
       if (.not.inp_i(solve_maxsub)) solve_maxsub = 10
       if (.not. rtdb_put(rtdb,'gradients:solve_maxsub',mt_int,
     &           1,solve_maxsub))
     &  call errquit(pname//'failed to write solve_maxsub',0,RTDB_ERR)
c 
      else if (inp_compare(.false.,test,'solve_maxit')) then
       if (.not.inp_i(solve_maxit)) solve_maxit = 250
       if (.not. rtdb_put(rtdb,'gradients:solve_maxit',mt_int,
     &           1,solve_maxit))
     &  call errquit(pname//'failed to write solve_maxit',0,RTDB_ERR)
c
c     Print
      else if (inp_compare(.false.,test,'print')) then
        call util_print_input(rtdb,'gradients')
c
c     END
      else if (inp_compare(.false.,test,'end')) then
        goto 20
c
      else
        call errquit(pname//'unknown directive',0, INPUT_ERR)
      endif
      goto 10
c
c ------
c Return
c ------
c
 20   return
      end
