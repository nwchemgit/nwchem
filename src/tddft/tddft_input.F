      subroutine tddft_input(rtdb)
c
c $Id: tddft_input.F,v 1.1 2002-03-22 02:10:37 sohirata Exp $
c
c Input parser for TDDFT called by input/input_parse.F.
c Also sets default values for input parameters.
c Modified from mp2_grad/mp2_input.F by So Hirata Feb, 2002.
c
c     TDDFT
c        [(TDA||NOTDA) default NOTDA]
c        [NROOTS <integer nroots default 1>]
c        [MAXVECS <integer maxvecs default 1000>]
c        [(SINGLET||NOSINGLET) default SINGLET]
c        [(TRIPLET||NOTRIPLET) default TRIPLET]
c        [THRESH <double thresh default 1e-5>]
c        [MAXITER <integer maxiter default 100>]
c        [TARGET <integer target default 1>]
c        [ALGORITHM <integer algorithm default 1>]
c        [PRINT (none||low||medium||high||debug)
c          <string list_of_names ...>]
c     END
c     
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb
      character*20 test
      integer nroots,maxvecs,maxiter,target,algorithm
      integer ihuge
      parameter (ihuge=1000000000)
      logical singlet,triplet
      logical tda
      double precision thresh
c
c ------------------
c Set default values
c ------------------
c
c     TDA (Tamm-Dancoff approximation)
c
      tda=.false.
      if (.not.rtdb_put(rtdb,'tddft:tda',mt_log,1,tda))
     1  call errquit('tddft_input: failed writing to rtdb',0)
c
c     NROOTS (the number of excited states sought) 
c
      nroots=1
      if (.not.rtdb_put(rtdb,'tddft:nroots',mt_int,1,nroots))
     1  call errquit('tddft_input: failed writing to rtdb',0)
c
c     MAXVECS (the maximum number of trial vectors)
c
      maxvecs=ihuge
      if (.not.rtdb_put(rtdb,'tddft:maxvecs',mt_int,1,maxvecs))
     1  call errquit('tddft_input: failed writing to rtdb',0)
c
c     SINGLET (singlet excited states sought for a closed shell)
c
      singlet=.true.
      if (.not.rtdb_put(rtdb,'tddft:singlet',mt_log,1,singlet))
     1  call errquit('tddft_input: failed writing to rtdb',0)
c
c     TRIPLET (triplet excited states sought for a closed shell)
c
      triplet=.true.
      if (.not.rtdb_put(rtdb,'tddft:triplet',mt_log,1,triplet))
     1  call errquit('tddft_input: failed writing to rtdb',0)
c
c     THRESH (convergence threshold for Davidson iteration)
c
      thresh=1.0d-5
      if (.not.rtdb_put(rtdb,'tddft:thresh',mt_dbl,1,thresh))
     1  call errquit('tddft_input: failed writing to rtdb',0)
c
c     MAXITER (the maximum number of Davidson iterations)
c
      maxiter=100
      if (.not.rtdb_put(rtdb,'tddft:maxiter',mt_int,1,maxiter))
     1  call errquit('tddft_input: failed writing to rtdb',0)
c
c     TARGET (the target excited state for, e.g., geometry optimization)
c
      target=1
      if (.not.rtdb_put(rtdb,'tddft:target',mt_int,1,target))
     1  call errquit('tddft_input: failed writing to rtdb',0)
c
c     ALGORITHM (algorithm: =1 multiple fock build; =2: conventional)
c
      algorithm=1
      if (.not.rtdb_put(rtdb,'tddft:algorithm',mt_int,1,algorithm))
     1  call errquit('tddft_input: failed writing to rtdb',0)
c
c ----------
c Read input
c ----------
c
 10   if (.not. inp_read()) 
     1  call errquit('tddft_input: failed reading input',0)
      if (.not. inp_a(test)) 
     1  call errquit('tddft_input: failed reading keyword',0)
c
c     TDA (Tamm-Dancoff approximation)
c
      if (inp_compare(.false.,test,'tda')) then
        tda=.true.
        if (.not.rtdb_put(rtdb,'tddft:tda',mt_log,1,tda))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     NOTDA (Full linear response)
c
      else if (inp_compare(.false.,test,'notda')) then
        tda=.false.
        if (.not.rtdb_put(rtdb,'tddft:tda',mt_log,1,tda))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     NROOTS (the number of excited states sought)
c
      else if (inp_compare(.false.,test,'nroots')) then
        if (.not.inp_i(nroots)) then
          write(LuOut,*) 'tddft_input: nroots value not found; ',
     1      'default value of 1 will be used'
          nroots=1
        endif
        if (.not.rtdb_put(rtdb,'tddft:nroots',mt_int,1,nroots))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     MAXVECS (the maximum number of trial vectors)
c
      else if (inp_compare(.false.,test,'maxvecs')) then
        if (.not.inp_i(maxvecs)) then
          write(LuOut,*) 'tddft_input: maxvecs value not found; ',
     1      'will be determined based on available memory size'
          maxvecs=ihuge
        endif
        if (.not.rtdb_put(rtdb,'tddft:maxvecs',mt_int,1,maxvecs))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     NOSINGLET (no singlet excited states sought for a closed shell)
c
      else if (inp_compare(.false.,test,'nosinglet')) then
        singlet=.false.
        if (.not.rtdb_put(rtdb,'tddft:singlet',mt_log,1,singlet))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     SINGLET (singlet excited states sought for a closed shell)
c
      else if (inp_compare(.false.,test,'singlet')) then
        singlet=.true.
        if (.not.rtdb_put(rtdb,'tddft:singlet',mt_log,1,singlet))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     NOTRIPLET (no triplet excited states sought for a closed shell)
c
      else if (inp_compare(.false.,test,'notriplet')) then
        triplet=.false.
        if (.not.rtdb_put(rtdb,'tddft:triplet',mt_log,1,triplet))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     TRIPLET (triplet excited states sought for a closed shell)
c
      else if (inp_compare(.false.,test,'triplet')) then
        triplet=.true.
        if (.not.rtdb_put(rtdb,'tddft:triplet',mt_log,1,triplet))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     THRESH (convergence threshold for Davidson iteration)
c
      else if (inp_compare(.false.,test,'thresh')) then
        if (.not.inp_f(thresh)) then
          write(LuOut,*) 'tddft_input: thresh value not found; ',
     1      'default value of 1e-5 will be used'
          thresh=1.0d-5
        endif
        if (.not.rtdb_put(rtdb,'tddft:thresh',mt_dbl,1,thresh))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     MAXITER (the maximum number of Davidson iterations)
c
      else if (inp_compare(.false.,test,'maxiter')) then
        if (.not.inp_i(maxiter)) then
          write(LuOut,*) 'tddft_input: maxiter value not found; ',
     1      'default value of 100 will be used'
          maxiter=100
        endif
        if (.not.rtdb_put(rtdb,'tddft:maxiter',mt_int,1,maxiter))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     TARGET (the target excited state for, e.g., geometry optimization)
c
      else if (inp_compare(.false.,test,'target')) then
        if (.not.inp_i(target)) then
          write(LuOut,*) 'tddft_input: target value not found; ',
     1      'default value of 1 will be used'
          target=1
        endif
        if (target.gt.nroots) call errquit
     1    ('tddft_input: an illegal value for target',target)
        if (.not.rtdb_put(rtdb,'tddft:target',mt_int,1,target))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     ALGORITHM (algorithm: =1 multiple fock build; =2: conventional)
c
      else if (inp_compare(.false.,test,'algorithm')) then
        if (.not.inp_i(algorithm)) then
          write(LuOut,*) 'tddft_input: algorithm value not found; ',
     1      'default value of 1 will be used'
          algorithm=1
        endif
        if (.not.rtdb_put(rtdb,'tddft:algorithm',mt_int,1,algorithm))
     1    call errquit('tddft_input: failed writing to rtdb',0)
c
c     PRINT
c
      else if (inp_compare(.false.,test,'print')) then
        call util_print_input(rtdb,'tddft')
c
c     END
c
      else if (inp_compare(.false.,test,'end')) then
        goto 20
      else
        call errquit('tddft_input: unknown directive',0)
      endif
      goto 10
c
c ------
c Return
c ------
c
 20   return
      end
